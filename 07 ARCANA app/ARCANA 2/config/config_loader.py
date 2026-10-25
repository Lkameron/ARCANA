"""
ARCANA config loader — reads project.yaml and exposes helpers used by
the backend and the AI agent.
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml")

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "project.yaml"

_cache: dict | None = None


def load(path: Path = CONFIG_PATH) -> dict:
    global _cache
    if _cache is None:
        with open(path) as f:
            _cache = yaml.safe_load(f)
    return _cache


def reload(path: Path = CONFIG_PATH) -> dict:
    global _cache
    _cache = None
    return load(path)


# ── Convenience accessors ─────────────────────────────────────────────────────

def project_name() -> str:
    return load().get("project", {}).get("name", "ARCANA Project")

def dataset_id() -> str:
    return load().get("project", {}).get("dataset_id", "")

def description() -> str:
    return load().get("project", {}).get("description", "")

def species() -> str:
    return load().get("species", {}).get("name", "human")

def genome() -> str:
    return load().get("species", {}).get("genome", "GRCh38")

def comparisons() -> list[dict]:
    return load().get("comparisons", [])

def comparison_names() -> list[str]:
    return [c["name"] for c in comparisons()]

def first_comparison() -> str:
    names = comparison_names()
    return names[0] if names else ""

def traits() -> dict:
    return load().get("traits", {"continuous": [], "binary": []})

def condition_col() -> str:
    return load().get("metadata", {}).get("condition_col", "condition")

def sample_id_col() -> str:
    return load().get("metadata", {}).get("sample_id_col", "sample_id")

def results_dir() -> Path:
    return ROOT / load().get("output", {}).get("results_dir", "results")

def counts_path() -> Path:
    return ROOT / load().get("data", {}).get("counts_matrix", "data/processed/counts_matrix.tsv")

def metadata_path() -> Path:
    return ROOT / load().get("data", {}).get("metadata", "data/metadata/sample_metadata.csv")

def de_thresholds() -> dict:
    return load().get("de_thresholds", {"padj": 0.05, "log2fc": 1.0})

def wgcna_settings() -> dict:
    return load().get("wgcna", {})

def lncrna_biotypes() -> list[str]:
    return load().get("lncrna", {}).get("biotypes", ["lncRNA"])

def lncrna_enabled() -> bool:
    return load().get("lncrna", {}).get("enabled", True)


# ── Environment variable dict (passed to R / Python scripts) ──────────────────

def as_env() -> dict[str, str]:
    """
    Returns a dict of environment variables that every pipeline script reads.
    Passed via subprocess env so no script hard-codes anything.
    """
    cfg = load()
    t = traits()
    wgcna = wgcna_settings()
    thresh = de_thresholds()
    comps = comparisons()

    return {
        "ARCANA_PROJECT_NAME":      project_name(),
        "ARCANA_DATASET_ID":        dataset_id(),
        "ARCANA_DESCRIPTION":       description().strip(),
        "ARCANA_SPECIES":           species(),
        "ARCANA_GENOME":            genome(),
        "ARCANA_ENSEMBL_VERSION":   str(cfg.get("species", {}).get("ensembl_version", 109)),
        "ARCANA_COUNTS_MATRIX":     str(counts_path()),
        "ARCANA_METADATA":          str(metadata_path()),
        "ARCANA_SAMPLE_ID_COL":     sample_id_col(),
        "ARCANA_CONDITION_COL":     condition_col(),
        "ARCANA_COMPARISONS":       ",".join(c["name"] for c in comps),
        "ARCANA_COMPARISONS_JSON":  json.dumps(comps),
        "ARCANA_TRAITS_CONTINUOUS": ",".join(t.get("continuous") or []),
        "ARCANA_TRAITS_BINARY":     ",".join(t.get("binary") or []),
        "ARCANA_RESULTS_DIR":       str(results_dir()),
        "ARCANA_PADJ":              str(thresh.get("padj", 0.05)),
        "ARCANA_LOG2FC":            str(thresh.get("log2fc", 1.0)),
        "ARCANA_WGCNA_SOFT_POWER":  str(wgcna.get("soft_power") or ""),
        "ARCANA_WGCNA_MIN_MODULE":  str(wgcna.get("min_module_size", 30)),
        "ARCANA_WGCNA_MERGE_CUT":   str(wgcna.get("merge_cut_height", 0.25)),
        "ARCANA_WGCNA_N_FOCUS":     str(wgcna.get("n_focus_modules", 5)),
        "ARCANA_LNCRNA_ENABLED":    str(lncrna_enabled()).upper(),
        "ARCANA_LNCRNA_BIOTYPES":   ",".join(lncrna_biotypes()),
        "ARCANA_DE_COVARIATES":     ",".join(cfg.get("de_covariates") or []),
        "ARCANA_FIGURES_DPI":       str(cfg.get("output", {}).get("figures_dpi", 300)),
    }


# ── AI system prompt builder ───────────────────────────────────────────────────

def build_ai_system_prompt() -> str:
    """
    Builds the AI agent system prompt dynamically from the project config,
    so the AI knows the researcher's specific dataset without any hardcoding.
    """
    cfg = load()
    comp_list = "\n".join(
        f"  - {c['name']}: {c.get('label', c['name'])}"
        for c in comparisons()
    )
    t = traits()
    cont_traits = ", ".join(t.get("continuous") or []) or "none defined"
    bin_traits  = ", ".join(t.get("binary") or []) or "none defined"

    return f"""You are ARCANA-AI, an expert bioinformatics assistant embedded in the ARCANA
lncRNA co-expression analysis platform. You have full access to this researcher's
analysis results and can read, summarise, and explain them.

## Current project
- Name: {project_name()}
- Dataset: {dataset_id()}
- Description: {description().strip()}
- Species: {species()} ({genome()})
- Condition column: {condition_col()}

## Comparisons analysed
{comp_list if comp_list.strip() else "  (none configured yet)"}

## Trait correlations
- Continuous traits: {cont_traits}
- Binary traits: {bin_traits}

## Analysis pipeline
The pipeline runs these steps in order:
1. QC & normalisation (edgeR/DESeq2 normalisation, low-count filtering)
2. Differential expression (DESeq2 with apeglm shrinkage)
3. lncRNA prioritisation (Ensembl biotype annotation, expression filter)
4. WGCNA co-expression network (soft thresholding, module detection)
5. Hub gene prioritisation (kME module membership + intramodular connectivity)
6. Enrichment analysis (GO Biological Process via clusterProfiler)
7. ML ranking (XGBoost with SHAP feature importance)
8. Publication figures
9. Automated HTML report

## Your role
- Answer questions about the results, methods, and biological interpretation
- Read result files when asked (use the list_files and read_table tools)
- Summarise DE results, hub genes, and enrichment findings
- Explain statistical concepts in plain language when asked
- Do NOT fabricate gene names, statistics, or pathway results — always read from files

## Important
- All analysis runs locally on the researcher's machine
- No data is sent externally
- If Ollama is not available, politely explain the researcher needs to install it
"""
