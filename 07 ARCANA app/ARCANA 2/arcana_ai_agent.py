"""
ARCANA-AI — Local bioinformatics agent for RNA-seq co-expression analysis.

Requires Ollama running locally:
  brew install ollama
  ollama serve &
  ollama pull llama3.2   (or mistral, gemma2, etc.)

Run in terminal:
  python3.11 scripts/arcana_ai_agent.py

Or import and use the ArcanaAgent class from your own code.
"""

import json
import subprocess
import sys
import os
import gzip
import shutil
from pathlib import Path
from typing import Any

import requests
import pandas as pd

# ── Config ────────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS      = PROJECT_ROOT / "results"
DATA         = PROJECT_ROOT / "data"

OLLAMA_URL   = os.environ.get("OLLAMA_URL",  "http://localhost:11434")
OLLAMA_MODEL = os.environ.get("ARCANA_MODEL", "llama3.2")

SYSTEM_PROMPT = """You are ARCANA-AI, a specialist bioinformatics assistant for RNA-seq analysis and weighted gene co-expression network analysis (WGCNA).

You run entirely locally on the user's machine. You have tools to:
• Inspect FASTQ / FASTA / BAM files and report their basic statistics
• Run FastQC quality control on raw sequencing files
• Execute ARCANA pipeline steps (normalisation, DESeq2, WGCNA, ML ranking)
• Read result tables (DE results, hub gene rankings, module-trait correlations)
• Summarise differential expression results for any comparison
• List files in the project directory

Current project: Type 2 Diabetes lncRNA co-expression analysis of pancreatic islet bulk RNA-seq (GSE164416).
  - 36 donors: T2D (n=12), IGT (n=9), non-diabetic ND (n=15)
  - 3 comparisons: T2D_vs_ND, IGT_vs_ND, T2D_vs_IGT
  - 30 WGCNA modules; 5 focus modules: brown, yellow, green, red, black
  - 604 significant DE genes (T2D_vs_ND, padj<0.05)
  - 794 hub lncRNAs identified across focus modules

Key result paths:
  results/de/tables/           — DESeq2 tables
  results/wgcna/hubgenes/tables/ — hub gene rankings
  results/wgcna/module_trait/  — module-trait correlations
  results/ml/                  — XGBoost / SHAP outputs
  data/raw/                    — raw FASTQ / BAM inputs

Rules:
- Always call the appropriate tool before making claims about specific numbers or genes.
- Explain biological significance of results in plain language.
- If a pipeline step has not been run yet, tell the user what command to run.
- Be concise. Use bullet points for lists of genes or results.
"""

# ── Tool registry ─────────────────────────────────────────────────────────────

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": "List files in a project directory. Use to explore FASTQ, BAM, CSV, or result files.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Relative path from project root, e.g. 'data/raw', 'results/de/tables'"
                    },
                    "pattern": {
                        "type": "string",
                        "description": "Optional glob pattern, e.g. '*.fastq.gz', '*.csv'"
                    }
                },
                "required": ["path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "read_table",
            "description": "Read a TSV or CSV results table and return the first N rows. Use for DE results, hub gene rankings, module-trait tables, or any tabular output.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Relative path to TSV/CSV file from project root"
                    },
                    "rows": {
                        "type": "integer",
                        "description": "Number of rows to return (default 20)"
                    },
                    "columns": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Specific column names to include (all if omitted)"
                    }
                },
                "required": ["path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "summarize_de_results",
            "description": "Summarise differential expression results for a specific comparison. Returns counts, top up/down-regulated genes, and lncRNA fraction.",
            "parameters": {
                "type": "object",
                "properties": {
                    "comparison": {
                        "type": "string",
                        "description": "Comparison name: T2D_vs_ND, IGT_vs_ND, or T2D_vs_IGT"
                    },
                    "padj_threshold": {
                        "type": "number",
                        "description": "Adjusted p-value threshold (default 0.05)"
                    },
                    "lncrna_only": {
                        "type": "boolean",
                        "description": "Only return lncRNA genes (default false)"
                    }
                },
                "required": ["comparison"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_hub_lncrnas",
            "description": "Get top hub lncRNA genes from WGCNA analysis, ranked by hub score. Can filter by module colour.",
            "parameters": {
                "type": "object",
                "properties": {
                    "module": {
                        "type": "string",
                        "description": "Module colour to filter (e.g. brown, yellow, green, red, black). Leave blank for all focus modules."
                    },
                    "n": {
                        "type": "integer",
                        "description": "Number of top genes to return (default 15)"
                    }
                }
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "file_info",
            "description": "Get information about a FASTQ, FASTA, or BAM file: file size, estimated read count, compression status.",
            "parameters": {
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the file (absolute or relative to project root)"
                    }
                },
                "required": ["file_path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_fastqc",
            "description": "Run FastQC quality control on one or more FASTQ files. FastQC must be installed (brew install fastqc).",
            "parameters": {
                "type": "object",
                "properties": {
                    "fastq_path": {
                        "type": "string",
                        "description": "Path to FASTQ or FASTQ.gz file (or directory containing them)"
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Output directory for FastQC report (default: results/qc/fastqc)"
                    }
                },
                "required": ["fastq_path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_pipeline_step",
            "description": "Run a specific ARCANA pipeline analysis step. Each step runs the corresponding script in scripts/.",
            "parameters": {
                "type": "object",
                "properties": {
                    "step": {
                        "type": "string",
                        "description": "Pipeline step identifier",
                        "enum": [
                            "01_qc",
                            "02_normalize",
                            "03_deseq2",
                            "04_wgcna",
                            "05_hub_genes",
                            "08_ml_models",
                            "10_report_writer"
                        ]
                    },
                    "args": {
                        "type": "string",
                        "description": "Optional extra arguments to pass to the script"
                    }
                },
                "required": ["step"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "module_trait_summary",
            "description": "Summarise which WGCNA modules are significantly correlated with clinical traits (T2D status, HbA1c, BMI).",
            "parameters": {
                "type": "object",
                "properties": {
                    "trait": {
                        "type": "string",
                        "description": "Specific trait to focus on (e.g. status_T2D, hba1c, bmi). Leave blank for all."
                    },
                    "r_threshold": {
                        "type": "number",
                        "description": "Minimum absolute Pearson r to report (default 0.4)"
                    }
                }
            }
        }
    }
]

# ── Tool implementations ───────────────────────────────────────────────────────

def _resolve(path: str) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return PROJECT_ROOT / path


def tool_list_files(path: str, pattern: str = "*") -> dict:
    p = _resolve(path)
    if not p.exists():
        return {"error": f"Path not found: {p}"}
    files = sorted(p.glob(pattern or "*"))
    entries = []
    for f in files[:60]:
        size = f.stat().st_size if f.is_file() else 0
        entries.append({
            "name": f.name,
            "type": "dir" if f.is_dir() else "file",
            "size_kb": round(size / 1024, 1) if size else None
        })
    return {"path": str(p), "count": len(files), "entries": entries[:60]}


def tool_read_table(path: str, rows: int = 20, columns: list = None) -> dict:
    p = _resolve(path)
    if not p.exists():
        return {"error": f"File not found: {p}"}
    try:
        sep = "\t" if p.suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(p, sep=sep, nrows=rows)
        if columns:
            df = df[[c for c in columns if c in df.columns]]
        return {
            "path": str(p),
            "shape": list(df.shape),
            "columns": list(df.columns),
            "rows": df.head(rows).to_dict(orient="records")
        }
    except Exception as e:
        return {"error": str(e)}


def tool_summarize_de(comparison: str, padj_threshold: float = 0.05, lncrna_only: bool = False) -> dict:
    suffix = "lncRNA_only" if lncrna_only else "all_genes"
    candidates = [
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}_{suffix}.tsv",
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}_{suffix}.csv",
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}.tsv",
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}.csv",
    ]
    path = next((c for c in candidates if c.exists()), None)
    if not path:
        return {"error": f"DE results not found for {comparison}. Run step 03_deseq2 first.",
                "tried": [str(c) for c in candidates]}
    try:
        sep = "\t" if path.suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(path, sep=sep)
        padj_col = next((c for c in ["padj","adj.P.Val","FDR"] if c in df.columns), None)
        lfc_col  = next((c for c in ["log2FoldChange","logFC","LFC"] if c in df.columns), None)
        sym_col  = next((c for c in ["gene_symbol","gene_name","gene_id","gene_symbol.x"] if c in df.columns), None)
        bio_col  = next((c for c in ["biotype","gene_biotype","biotype.x"] if c in df.columns), None)

        sig = df[df[padj_col] < padj_threshold] if padj_col else df
        up  = sig[sig[lfc_col] > 0] if lfc_col else pd.DataFrame()
        dn  = sig[sig[lfc_col] < 0] if lfc_col else pd.DataFrame()

        lnc_count = 0
        if bio_col is not None:
            lnc_count = sig[sig[bio_col].str.contains("lncRNA|lnc|antisense|lincRNA", case=False, na=False)].shape[0]

        top_up = up.nlargest(5, lfc_col)[sym_col].tolist() if (lfc_col and sym_col and not up.empty) else []
        top_dn = dn.nsmallest(5, lfc_col)[sym_col].tolist() if (lfc_col and sym_col and not dn.empty) else []

        return {
            "comparison": comparison,
            "total_genes_tested": len(df),
            "significant_padj_lt": padj_threshold,
            "n_significant": len(sig),
            "n_upregulated": len(up),
            "n_downregulated": len(dn),
            "n_lncrna_significant": lnc_count,
            "top_upregulated_genes": top_up,
            "top_downregulated_genes": top_dn,
        }
    except Exception as e:
        return {"error": str(e)}


def tool_get_hub_lncrnas(module: str = "", n: int = 15) -> dict:
    candidates = [
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_hub_genes_ranked_focus_modules.csv",
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_top20_hub_genes_per_focus_module.csv",
    ]
    path = next((c for c in candidates if c.exists()), None)
    if not path:
        return {"error": "Hub gene table not found. Run step 05_hub_genes first."}
    try:
        df = pd.read_csv(path)
        bio_col = next((c for c in ["biotype.x","biotype","gene_biotype"] if c in df.columns), None)
        sym_col = next((c for c in ["gene_symbol.x","gene_symbol","gene_id"] if c in df.columns), None)
        mod_col = next((c for c in ["module_color","module"] if c in df.columns), None)
        scr_col = next((c for c in ["hub_score_lncRNA","hub_score_all_genes","own_kME"] if c in df.columns), None)

        if bio_col:
            df = df[df[bio_col].str.contains("lncRNA|lnc|antisense|lincRNA", case=False, na=False)]
        if module and mod_col:
            df = df[df[mod_col].str.lower() == module.lower()]

        if scr_col:
            df = df.sort_values(scr_col, ascending=False)

        cols = [c for c in [sym_col, mod_col, scr_col, bio_col] if c]
        top = df[cols].head(n).to_dict(orient="records") if cols else df.head(n).to_dict(orient="records")

        return {
            "module_filter": module or "all focus modules",
            "n_lncrna_total": len(df),
            "top_hub_lncrnas": top
        }
    except Exception as e:
        return {"error": str(e)}


def tool_file_info(file_path: str) -> dict:
    p = _resolve(file_path)
    if not p.exists():
        return {"error": f"File not found: {p}"}

    info = {
        "path": str(p),
        "size_mb": round(p.stat().st_size / (1024**2), 2),
        "suffix": p.suffix,
    }

    suffix = p.suffix.lower()
    # FASTQ read count estimate
    if suffix in (".fastq", ".fq"):
        try:
            result = subprocess.run(["wc", "-l", str(p)], capture_output=True, text=True, timeout=10)
            lines = int(result.stdout.strip().split()[0])
            info["estimated_reads"] = lines // 4
            info["format"] = "FASTQ (uncompressed)"
        except Exception:
            pass
    elif suffix == ".gz":
        if ".fastq" in p.name or ".fq" in p.name:
            try:
                # Count first 40000 bytes to estimate
                with gzip.open(p, "rt") as f:
                    sample = sum(1 for _ in zip(range(40000), f))
                info["format"] = "FASTQ (gzipped)"
                info["note"] = f"Sampled {sample} lines from start; full read count requires full scan"
            except Exception:
                info["format"] = "FASTQ (gzipped)"
    elif suffix == ".bam":
        info["format"] = "BAM"
        if shutil.which("samtools"):
            try:
                result = subprocess.run(
                    ["samtools", "flagstat", str(p)],
                    capture_output=True, text=True, timeout=30
                )
                info["samtools_flagstat"] = result.stdout[:600]
            except Exception:
                info["note"] = "samtools available but flagstat failed"
        else:
            info["note"] = "Install samtools for BAM statistics: brew install samtools"
    elif suffix in (".fa", ".fasta", ".fna"):
        info["format"] = "FASTA"
        try:
            result = subprocess.run(
                ["grep", "-c", "^>", str(p)],
                capture_output=True, text=True, timeout=15
            )
            info["sequence_count"] = int(result.stdout.strip())
        except Exception:
            pass

    return info


def tool_run_fastqc(fastq_path: str, output_dir: str = "results/qc/fastqc") -> dict:
    if not shutil.which("fastqc"):
        return {
            "error": "FastQC not installed",
            "fix": "Run: brew install fastqc"
        }
    p = _resolve(fastq_path)
    out = _resolve(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if not p.exists():
        return {"error": f"File not found: {p}"}

    try:
        result = subprocess.run(
            ["fastqc", str(p), "-o", str(out), "--quiet"],
            capture_output=True, text=True, timeout=300
        )
        report = out / (p.stem.replace(".fastq","").replace(".fq","") + "_fastqc.html")
        return {
            "status": "success" if result.returncode == 0 else "error",
            "output_dir": str(out),
            "report": str(report) if report.exists() else "report not found",
            "stderr": result.stderr[-400:] if result.stderr else "",
        }
    except subprocess.TimeoutExpired:
        return {"error": "FastQC timed out after 5 minutes"}
    except Exception as e:
        return {"error": str(e)}


STEP_SCRIPTS = {
    "01_qc":          ("python3.11", "scripts/01_qc.py"),
    "02_normalize":   ("python3.11", "scripts/02_normalize.py"),
    "03_deseq2":      ("Rscript",    "scripts/03_deseq2.R"),
    "04_wgcna":       ("Rscript",    "scripts/04_wgcna.R"),
    "05_hub_genes":   ("Rscript",    "scripts/05_hub_genes.R"),
    "08_ml_models":   ("python3.11", "scripts/08_ml_models.py"),
    "10_report_writer":("python3.11","scripts/10_report_writer.py"),
}

def tool_run_pipeline_step(step: str, args: str = "") -> dict:
    if step not in STEP_SCRIPTS:
        return {"error": f"Unknown step: {step}. Valid: {list(STEP_SCRIPTS)}"}

    interp, script = STEP_SCRIPTS[step]
    script_path = PROJECT_ROOT / script
    if not script_path.exists():
        return {"error": f"Script not found: {script_path}"}

    exe = shutil.which(interp)
    if not exe:
        return {"error": f"{interp} not found in PATH"}

    cmd = [exe, str(script_path)]
    if args:
        cmd += args.split()

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=1800, cwd=str(PROJECT_ROOT)
        )
        tail = lambda s: s[-800:] if len(s) > 800 else s
        return {
            "step": step,
            "script": script,
            "returncode": result.returncode,
            "stdout_tail": tail(result.stdout),
            "stderr_tail": tail(result.stderr),
            "status": "success" if result.returncode == 0 else "failed"
        }
    except subprocess.TimeoutExpired:
        return {"error": f"Step {step} timed out after 30 minutes"}
    except Exception as e:
        return {"error": str(e)}


def tool_module_trait_summary(trait: str = "", r_threshold: float = 0.4) -> dict:
    corr_path = RESULTS / "wgcna" / "module_trait" / "wgcna_module_trait_correlations.csv"
    pval_path = RESULTS / "wgcna" / "module_trait" / "wgcna_module_trait_pvalues.csv"

    if not corr_path.exists():
        return {"error": "Module-trait correlations not found. Run step 04_wgcna first."}

    try:
        corr = pd.read_csv(corr_path, index_col=0)
        pval = pd.read_csv(pval_path, index_col=0) if pval_path.exists() else None

        results = []
        cols = [c for c in corr.columns if (not trait or trait.lower() in c.lower())]
        for col in cols:
            for module in corr.index:
                r = corr.loc[module, col]
                if abs(r) >= r_threshold:
                    p = pval.loc[module, col] if pval is not None and col in pval.columns else None
                    results.append({
                        "module": module,
                        "trait": col,
                        "r": round(float(r), 3),
                        "padj": round(float(p), 4) if p is not None else None,
                        "direction": "positive" if r > 0 else "negative"
                    })

        results.sort(key=lambda x: abs(x["r"]), reverse=True)
        return {
            "r_threshold": r_threshold,
            "trait_filter": trait or "all",
            "n_significant_associations": len(results),
            "associations": results[:30]
        }
    except Exception as e:
        return {"error": str(e)}


# ── Tool dispatcher ───────────────────────────────────────────────────────────

def dispatch_tool(name: str, arguments: dict) -> str:
    try:
        if name == "list_files":
            result = tool_list_files(**arguments)
        elif name == "read_table":
            result = tool_read_table(**arguments)
        elif name == "summarize_de_results":
            result = tool_summarize_de(**arguments)
        elif name == "get_hub_lncrnas":
            result = tool_get_hub_lncrnas(**arguments)
        elif name == "file_info":
            result = tool_file_info(**arguments)
        elif name == "run_fastqc":
            result = tool_run_fastqc(**arguments)
        elif name == "run_pipeline_step":
            result = tool_run_pipeline_step(**arguments)
        elif name == "module_trait_summary":
            result = tool_module_trait_summary(**arguments)
        else:
            result = {"error": f"Unknown tool: {name}"}
    except Exception as e:
        result = {"error": f"Tool execution error: {e}"}

    return json.dumps(result, default=str)


# ── Ollama agent loop ─────────────────────────────────────────────────────────

def check_ollama() -> bool:
    try:
        r = requests.get(f"{OLLAMA_URL}/api/tags", timeout=3)
        return r.status_code == 200
    except Exception:
        return False


def get_available_models() -> list:
    try:
        r = requests.get(f"{OLLAMA_URL}/api/tags", timeout=3)
        if r.ok:
            return [m["name"] for m in r.json().get("models", [])]
    except Exception:
        pass
    return []


def chat_turn(messages: list, stream_callback=None) -> str:
    """
    Send messages to Ollama, handle tool calls, return final text response.
    stream_callback(chunk: str) is called as text streams in (optional).
    """
    while True:
        payload = {
            "model": OLLAMA_MODEL,
            "messages": messages,
            "tools": TOOLS,
            "stream": False,
        }

        try:
            resp = requests.post(
                f"{OLLAMA_URL}/api/chat",
                json=payload,
                timeout=120
            )
            resp.raise_for_status()
        except requests.exceptions.ConnectionError:
            return "❌ Cannot reach Ollama. Is it running? Try: ollama serve"
        except Exception as e:
            return f"❌ Ollama error: {e}"

        data = resp.json()
        msg  = data.get("message", {})
        tool_calls = msg.get("tool_calls", [])

        # No tool calls → final text response
        if not tool_calls:
            content = msg.get("content", "")
            if stream_callback and content:
                stream_callback(content)
            messages.append({"role": "assistant", "content": content})
            return content

        # Has tool calls → execute each and loop
        messages.append({"role": "assistant", "content": msg.get("content",""), "tool_calls": tool_calls})

        for tc in tool_calls:
            fn   = tc.get("function", {})
            name = fn.get("name", "")
            args = fn.get("arguments", {})
            if isinstance(args, str):
                try:
                    args = json.loads(args)
                except Exception:
                    args = {}

            if stream_callback:
                stream_callback(f"\n⚙️  *{name}({', '.join(f'{k}={v!r}' for k,v in args.items())})*\n")

            result_str = dispatch_tool(name, args)

            messages.append({
                "role": "tool",
                "content": result_str,
            })


# ── Terminal chat interface ───────────────────────────────────────────────────

def run_terminal():
    global OLLAMA_MODEL
    print("\n" + "═"*60)
    print("  ARCANA-AI  —  Local RNA-seq Bioinformatics Agent")
    print("═"*60)

    if not check_ollama():
        print("\n❌  Ollama is not running.")
        print("    Start it with:  ollama serve")
        print(f"    Then pull a model:  ollama pull {OLLAMA_MODEL}")
        sys.exit(1)

    models = get_available_models()
    if OLLAMA_MODEL not in models:
        if models:
            OLLAMA_MODEL = models[0]
            print(f"⚠️  Model '{OLLAMA_MODEL}' not found. Using: {OLLAMA_MODEL}")
        else:
            print(f"\n❌  No models pulled yet.")
            print(f"    Run:  ollama pull llama3.2")
            sys.exit(1)

    print(f"\n✓  Connected to Ollama  ·  Model: {OLLAMA_MODEL}")
    print(f"   Project: {PROJECT_ROOT.name}")
    print(f"\n   Ask me about your RNA-seq data, DE results, hub lncRNAs,")
    print(f"   WGCNA modules, or drop FASTQ/BAM paths to analyse files.\n")
    print("   Type  exit  or  quit  to leave.\n")
    print("─"*60)

    messages = [{"role": "system", "content": SYSTEM_PROMPT}]

    while True:
        try:
            user_input = input("\nYou: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\nGoodbye.")
            break

        if not user_input:
            continue
        if user_input.lower() in ("exit", "quit", "q"):
            print("Goodbye.")
            break

        messages.append({"role": "user", "content": user_input})

        print("\nARCANA-AI: ", end="", flush=True)

        def print_chunk(text):
            print(text, end="", flush=True)

        response = chat_turn(messages, stream_callback=print_chunk)
        print()  # newline after response


# ── Programmatic API (used by FastAPI backend) ────────────────────────────────

class ArcanaAgent:
    """
    Stateless wrapper — create one instance, call .respond() per turn.
    Pass history in, get updated history + text response out.
    """

    def __init__(self, model: str = None):
        self.model = model or OLLAMA_MODEL
        self._prev_model = OLLAMA_MODEL

    def __enter__(self):
        global OLLAMA_MODEL
        self._prev_model = OLLAMA_MODEL
        OLLAMA_MODEL = self.model
        return self

    def __exit__(self, *_):
        global OLLAMA_MODEL
        OLLAMA_MODEL = self._prev_model

    def respond(self, history: list, user_message: str) -> tuple[str, list]:
        """
        history: list of {role, content} dicts (NOT including the new user message).
        Returns (assistant_text, updated_history).
        """
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        messages.extend(history)
        messages.append({"role": "user", "content": user_message})

        with self:
            response = chat_turn(messages)

        updated = history + [
            {"role": "user",      "content": user_message},
            {"role": "assistant", "content": response},
        ]
        return response, updated

    @staticmethod
    def is_available() -> bool:
        return check_ollama()

    @staticmethod
    def models() -> list:
        return get_available_models()


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_terminal()
