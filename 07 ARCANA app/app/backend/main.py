#!/usr/bin/env python3
"""
ARCANA — FastAPI backend
Serves the ARCANA UI and exposes REST + WebSocket endpoints that
drive the R/Python analysis pipeline.

Start:  uvicorn app.backend.main:app --reload --port 8000
"""

import asyncio
import json
import os
import sqlite3
import subprocess
import sys
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# ─── Paths ────────────────────────────────────────────────────────────────────

ROOT      = Path(__file__).resolve().parents[2]   # project root
SCRIPTS   = ROOT / "scripts"
RESULTS   = ROOT / "results"
DATA      = ROOT / "data"
FRONTEND  = Path(__file__).parent.parent / "frontend"
DB_PATH   = ROOT / "arcana_runs.db"

# ─── App ──────────────────────────────────────────────────────────────────────

app = FastAPI(title="ARCANA", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve the ARCANA UI at /
app.mount("/ui", StaticFiles(directory=str(FRONTEND), html=True), name="frontend")


@app.get("/")
def root():
    return FileResponse(str(FRONTEND / "index.html"))


# ─── Database ─────────────────────────────────────────────────────────────────

def get_db() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    conn = get_db()
    conn.execute("""
        CREATE TABLE IF NOT EXISTS runs (
            id          TEXT PRIMARY KEY,
            name        TEXT,
            status      TEXT DEFAULT 'queued',
            started_at  TEXT,
            finished_at TEXT,
            params      TEXT,
            steps_done  TEXT DEFAULT '[]',
            error       TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS run_files (
            id       INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id   TEXT REFERENCES runs(id),
            kind     TEXT,
            path     TEXT,
            label    TEXT
        )
    """)
    conn.commit()
    conn.close()


init_db()

# ─── WebSocket connection manager ─────────────────────────────────────────────

class ConnectionManager:
    def __init__(self):
        self.active: dict[str, list[WebSocket]] = {}

    async def connect(self, run_id: str, ws: WebSocket):
        await ws.accept()
        self.active.setdefault(run_id, []).append(ws)

    def disconnect(self, run_id: str, ws: WebSocket):
        self.active.get(run_id, []).remove(ws)

    async def broadcast(self, run_id: str, msg: dict):
        for ws in list(self.active.get(run_id, [])):
            try:
                await ws.send_json(msg)
            except Exception:
                pass


manager = ConnectionManager()


@app.websocket("/ws/run/{run_id}")
async def ws_run(run_id: str, websocket: WebSocket):
    await manager.connect(run_id, websocket)
    try:
        while True:
            await websocket.receive_text()   # keep alive
    except WebSocketDisconnect:
        manager.disconnect(run_id, websocket)


# ─── Pipeline step definitions ────────────────────────────────────────────────

PIPELINE_STEPS = [
    {"id": "qc",          "label": "QC + Normalise",          "script": "02_qc_filter_normalise.R",   "lang": "R"},
    {"id": "de",          "label": "Differential Expression",  "script": "03_differential_expression.R","lang": "R"},   # root
    {"id": "lncrna",      "label": "lncRNA Prioritisation",    "script": "04_lncRNA_prioritisation.R", "lang": "R"},
    {"id": "wgcna_prep",  "label": "WGCNA Inputs",             "script": "05A_prepare_wgcna_inputs.R", "lang": "R"},
    {"id": "wgcna_net",   "label": "WGCNA Network",            "script": "05B_wgcna_construct_network.R","lang": "R"},   # root
    # note: 05B is also in scripts/ as 03_wcgna_construct_network.R — kept flexible
    {"id": "hub",         "label": "Hub Gene Prioritisation",  "script": "05_hub_gene_prioritisation.R","lang": "R"},
    {"id": "enrichment",  "label": "Enrichment Analysis",      "script": "07_enrichment_analysis.R",   "lang": "R"},
    {"id": "ml",          "label": "ML Ranking + SHAP",        "script": "08_ml_models.py",             "lang": "python"},
    {"id": "figures",     "label": "Publication Figures",      "script": "09_figures_tables.R",         "lang": "R"},
    {"id": "report",      "label": "Report Writer",            "script": "10_report_writer.py",         "lang": "python"},
]

STEP_IDS = [s["id"] for s in PIPELINE_STEPS]


async def _run_script(step: dict, run_id: str, env: dict) -> bool:
    """Run a single pipeline script, streaming stdout to WebSocket."""
    script_path = SCRIPTS / step["script"]

    if not script_path.exists():
        # Try root directory fallback
        script_path = ROOT / step["script"]

    if not script_path.exists():
        await manager.broadcast(run_id, {
            "type": "log",
            "step": step["id"],
            "text": f"[SKIP] Script not found: {step['script']}",
        })
        return True   # treat missing optional scripts as success

    if step["lang"] == "R":
        cmd = ["Rscript", str(script_path)]
    else:
        cmd = [sys.executable, str(script_path)]

    await manager.broadcast(run_id, {
        "type": "step_start",
        "step": step["id"],
        "label": step["label"],
    })

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        cwd=str(ROOT),
        env={**os.environ, **env},
    )

    log_lines = []
    async for raw in proc.stdout:
        line = raw.decode(errors="replace").rstrip()
        log_lines.append(line)
        await manager.broadcast(run_id, {"type": "log", "step": step["id"], "text": line})

    await proc.wait()
    success = proc.returncode == 0

    await manager.broadcast(run_id, {
        "type": "step_done" if success else "step_error",
        "step": step["id"],
        "returncode": proc.returncode,
    })

    # Persist log
    log_file = RESULTS / "logs" / f"{run_id}_{step['id']}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("\n".join(log_lines))

    return success


async def _execute_run(run_id: str, steps: list[str], params: dict):
    """Background task: execute pipeline steps in order."""
    conn = get_db()
    conn.execute(
        "UPDATE runs SET status='running', started_at=? WHERE id=?",
        (datetime.utcnow().isoformat(), run_id),
    )
    conn.commit()
    conn.close()

    env = {
        "ARCANA_RUN_ID": run_id,
        "ARCANA_PARAMS": json.dumps(params),
    }

    steps_done = []
    for step in PIPELINE_STEPS:
        if step["id"] not in steps:
            continue
        ok = await _run_script(step, run_id, env)
        if ok:
            steps_done.append(step["id"])
            conn = get_db()
            conn.execute(
                "UPDATE runs SET steps_done=? WHERE id=?",
                (json.dumps(steps_done), run_id),
            )
            conn.commit()
            conn.close()
        else:
            conn = get_db()
            conn.execute(
                "UPDATE runs SET status='failed', finished_at=?, error=? WHERE id=?",
                (datetime.utcnow().isoformat(), f"Step {step['id']} failed", run_id),
            )
            conn.commit()
            conn.close()
            await manager.broadcast(run_id, {"type": "run_failed", "step": step["id"]})
            return

    conn = get_db()
    conn.execute(
        "UPDATE runs SET status='complete', finished_at=? WHERE id=?",
        (datetime.utcnow().isoformat(), run_id),
    )
    conn.commit()
    conn.close()

    await manager.broadcast(run_id, {"type": "run_complete", "run_id": run_id})

    # Index outputs
    _index_outputs(run_id)


def _index_outputs(run_id: str):
    """Record output files produced by this run into run_files table."""
    conn = get_db()
    for subdir, kind in [("figures/png", "figure"), ("figures/pdf", "figure_pdf"),
                          ("ml/plots", "ml_plot"), ("ml/tables", "ml_table"),
                          ("reports", "report")]:
        d = RESULTS / subdir
        if not d.exists():
            continue
        for f in d.iterdir():
            if f.is_file():
                conn.execute(
                    "INSERT OR IGNORE INTO run_files (run_id, kind, path, label) VALUES (?,?,?,?)",
                    (run_id, kind, str(f), f.name),
                )
    conn.commit()
    conn.close()


# ─── REST endpoints ───────────────────────────────────────────────────────────

class RunRequest(BaseModel):
    name: str = "Untitled run"
    steps: list[str] = STEP_IDS
    params: dict[str, Any] = {}


@app.post("/api/run")
async def create_run(req: RunRequest):
    run_id = str(uuid.uuid4())[:8]
    conn = get_db()
    conn.execute(
        "INSERT INTO runs (id, name, params) VALUES (?,?,?)",
        (run_id, req.name, json.dumps(req.params)),
    )
    conn.commit()
    conn.close()

    # Fire off background task
    asyncio.create_task(_execute_run(run_id, req.steps, req.params))

    return {"run_id": run_id, "status": "queued"}


@app.get("/api/runs")
def list_runs():
    conn = get_db()
    rows = conn.execute(
        "SELECT id, name, status, started_at, finished_at, steps_done FROM runs ORDER BY rowid DESC LIMIT 50"
    ).fetchall()
    conn.close()
    return [dict(r) for r in rows]


@app.get("/api/run/{run_id}")
def get_run(run_id: str):
    conn = get_db()
    row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if not row:
        raise HTTPException(404, "Run not found")
    files = conn.execute(
        "SELECT kind, path, label FROM run_files WHERE run_id=?", (run_id,)
    ).fetchall()
    conn.close()
    return {**dict(row), "files": [dict(f) for f in files]}


@app.get("/api/run/{run_id}/file")
def download_file(run_id: str, path: str):
    p = Path(path)
    if not p.exists():
        raise HTTPException(404, "File not found")
    return FileResponse(str(p), filename=p.name)


@app.get("/api/steps")
def list_steps():
    return PIPELINE_STEPS


# ─── Data endpoints ───────────────────────────────────────────────────────────

def _read_tsv_or_csv(path: Path, max_rows: int = 500) -> list:
    if not path.exists():
        return []
    sep = "\t" if path.suffix == ".tsv" else ","
    import csv
    return list(csv.DictReader(path.open(), delimiter=sep))[:max_rows]


def _first_existing(paths: list) -> Path | None:
    return next((p for p in paths if p.exists()), None)


@app.get("/api/results/de")
def get_de_results(comparison: str = "T2D_vs_ND", lncrna_only: bool = False):
    suffix = "lncRNA_only" if lncrna_only else "all_genes"
    candidates = [
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}_{suffix}.tsv",
        RESULTS / "de" / "tables" / f"DESeq2_{comparison}_significant.tsv",
        RESULTS / "de" / "tables" / "DESeq2_T2D_vs_ND_all_genes.tsv",
    ]
    p = _first_existing(candidates)
    return _read_tsv_or_csv(p) if p else []


@app.get("/api/results/de/comparisons")
def list_de_comparisons():
    tables_dir = RESULTS / "de" / "tables"
    if not tables_dir.exists():
        return []
    seen = set()
    out = []
    for f in sorted(tables_dir.iterdir()):
        if f.suffix == ".tsv" and "all_genes" in f.name:
            comp = f.stem.replace("DESeq2_", "").replace("_all_genes", "")
            if comp not in seen:
                seen.add(comp)
                out.append(comp)
    return out


@app.get("/api/results/hub_genes")
def get_hub_genes(lncrna_only: bool = False, module: str = "", limit: int = 200):
    p = _first_existing([
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_hub_genes_ranked_focus_modules.csv",
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_top20_hub_genes_per_focus_module.csv",
    ])
    rows = _read_tsv_or_csv(p, max_rows=5000) if p else []
    if lncrna_only:
        rows = [r for r in rows if "lncRNA" in r.get("biotype.x","") or "lncRNA" in r.get("biotype","")]
    if module:
        rows = [r for r in rows if r.get("module_color","") == module]
    return rows[:limit]


@app.get("/api/results/focus_modules")
def get_focus_modules():
    p = RESULTS / "wgcna" / "hubgenes" / "tables" / "focus_modules_selected_for_hub_analysis.csv"
    return _read_tsv_or_csv(p) if p.exists() else []


@app.get("/api/results/modules")
def get_modules():
    p = RESULTS / "wgcna" / "modules" / "wgcna_gene_module_assignments.csv"
    return _read_tsv_or_csv(p, max_rows=3000)


@app.get("/api/results/module_trait")
def get_module_trait():
    corr = _first_existing([RESULTS / "wgcna" / "module_trait" / "wgcna_module_trait_correlations.csv"])
    pval = _first_existing([RESULTS / "wgcna" / "module_trait" / "wgcna_module_trait_pvalues.csv"])
    assoc = _first_existing([RESULTS / "wgcna" / "module_trait" / "wgcna_module_status_associations_summary.csv"])
    return {
        "correlations": _read_tsv_or_csv(corr) if corr else [],
        "pvalues":      _read_tsv_or_csv(pval) if pval else [],
        "associations": _read_tsv_or_csv(assoc) if assoc else [],
    }


@app.get("/api/results/ml")
def get_ml_metrics():
    p = _first_existing([
        RESULTS / "ml" / "tables" / "binary_cv_metrics.csv",
        RESULTS / "ml" / "tables" / "binary_cv_metrics.tsv",
    ])
    return _read_tsv_or_csv(p) if p else []


@app.get("/api/results/enrichment")
def get_enrichment():
    p = _first_existing([
        RESULTS / "enrichment" / "GO_BP_enrichment_results.csv",
        RESULTS / "enrichment" / "enrichment_results.csv",
        RESULTS / "enrichment" / "tables" / "GO_BP_enrichment.csv",
    ])
    return _read_tsv_or_csv(p, max_rows=200) if p else []


@app.get("/api/results/summary")
def get_results_summary():
    """Single endpoint the UI calls on load to get all key numbers."""
    modules_p = RESULTS / "wgcna" / "modules" / "wgcna_module_sizes.csv"
    hub_p     = _first_existing([
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_hub_genes_ranked_focus_modules.csv",
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_top20_hub_genes_per_focus_module.csv",
    ])
    de_p      = _first_existing([
        RESULTS / "de" / "tables" / "DESeq2_T2D_vs_ND_significant.tsv",
    ])
    ml_p      = _first_existing([RESULTS / "ml" / "tables" / "binary_cv_metrics.csv"])

    module_rows = _read_tsv_or_csv(modules_p, max_rows=100) if modules_p.exists() else []
    hub_rows    = _read_tsv_or_csv(hub_p, max_rows=10) if hub_p else []
    de_sig      = _read_tsv_or_csv(de_p, max_rows=5000) if de_p else []
    ml_rows     = _read_tsv_or_csv(ml_p) if ml_p else []

    # top hub gene names (handle dot-suffixed column names from R merge)
    name_col = next((c for c in ["gene_symbol.x","gene_symbol","gene_name.x","gene_name","gene_id"] if hub_rows and c in hub_rows[0]), None)
    bio_col  = next((c for c in ["biotype.x","biotype","gene_biotype"] if hub_rows and c in hub_rows[0]), None)
    top_hubs = [
        {"symbol": r.get(name_col, r.get("gene_id","")), "biotype": r.get(bio_col,""), "module": r.get("module_color","")}
        for r in hub_rows[:10] if name_col
    ] if hub_rows else []

    # best ML model
    best_ml = {}
    if ml_rows:
        try:
            best_ml = max(ml_rows, key=lambda r: float(r.get("AUC_mean", 0) or 0))
        except Exception:
            pass

    return {
        "n_modules":        len(module_rows),
        "n_de_significant": len(de_sig),
        "top_hub_genes":    top_hubs,
        "best_ml_model":    best_ml.get("model", ""),
        "best_ml_auc":      best_ml.get("AUC_mean", ""),
        "comparisons_available": list_de_comparisons(),
        "figures_available": len(list_figures()),
    }


@app.get("/api/figures")
def list_figures():
    out = []
    search_dirs = [
        ("wgcna", RESULTS / "wgcna" / "plots"),
        ("wgcna", RESULTS / "wgcna" / "modules"),
        ("wgcna", RESULTS / "wgcna" / "module_trait"),
        ("de",    RESULTS / "de" / "plots"),
        ("ml",    RESULTS / "ml" / "plots"),
        ("qc",    RESULTS / "qc" / "plots"),
        ("figures", RESULTS / "figures"),
    ]
    for category, d in search_dirs:
        if d.exists():
            for f in sorted(d.iterdir()):
                if f.suffix.lower() in {".png", ".pdf", ".svg"}:
                    out.append({
                        "name":     f.name,
                        "path":     str(f),
                        "type":     f.suffix[1:],
                        "category": category,
                        "label":    f.stem.replace("_", " ").replace("wgcna ", "").title(),
                    })
    return out


@app.get("/api/figure")
def serve_figure(path: str):
    p = Path(path)
    if not p.exists():
        raise HTTPException(404)
    return FileResponse(str(p))


# ─── File upload ──────────────────────────────────────────────────────────────

@app.post("/api/upload/counts")
async def upload_counts(file: UploadFile = File(...)):
    dest = DATA / "processed" / "counts_matrix_upload.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(await file.read())
    return {"saved": str(dest), "filename": file.filename}


@app.post("/api/upload/metadata")
async def upload_metadata(file: UploadFile = File(...)):
    dest = DATA / "metadata" / "sample_metadata_upload.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(await file.read())
    return {"saved": str(dest), "filename": file.filename}


# ─── Document store (Phase 1b) ─────────────────────────────────────────────────

DOC_DB_PATH = ROOT / "arcana_docs.db"


def init_doc_db():
    conn = sqlite3.connect(str(DOC_DB_PATH))
    conn.execute("CREATE VIRTUAL TABLE IF NOT EXISTS docs USING fts5(path, title, body, tokenize='porter ascii')")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS doc_meta (
            path    TEXT PRIMARY KEY,
            title   TEXT,
            indexed TEXT
        )
    """)
    conn.commit()
    conn.close()


init_doc_db()


@app.post("/api/docs/index")
def index_documents(folder: str = str(ROOT / "dissertation")):
    """Index all PDFs/text files in folder into FTS5 store."""
    from pathlib import Path as P
    import importlib

    folder_path = P(folder)
    if not folder_path.exists():
        return {"indexed": 0, "error": "Folder not found"}

    conn = sqlite3.connect(str(DOC_DB_PATH))
    count = 0

    for f in folder_path.rglob("*"):
        if f.suffix.lower() == ".pdf":
            try:
                import pdfplumber
                with pdfplumber.open(f) as pdf:
                    body = "\n".join(p.extract_text() or "" for p in pdf.pages[:20])
                title = f.stem
            except Exception:
                continue
        elif f.suffix.lower() in {".md", ".txt", ".csv"}:
            body = f.read_text(errors="replace")[:50000]
            title = f.stem
        else:
            continue

        conn.execute("DELETE FROM docs WHERE path=?", (str(f),))
        conn.execute("INSERT INTO docs (path, title, body) VALUES (?,?,?)", (str(f), title, body))
        conn.execute(
            "INSERT OR REPLACE INTO doc_meta (path, title, indexed) VALUES (?,?,?)",
            (str(f), title, datetime.utcnow().isoformat()),
        )
        count += 1

    conn.commit()
    conn.close()
    return {"indexed": count}


@app.get("/api/docs/search")
def search_docs(q: str, limit: int = 10):
    conn = sqlite3.connect(str(DOC_DB_PATH))
    rows = conn.execute(
        "SELECT path, title, snippet(docs, 2, '<b>', '</b>', '...', 20) as excerpt "
        "FROM docs WHERE docs MATCH ? ORDER BY rank LIMIT ?",
        (q, limit),
    ).fetchall()
    conn.close()
    return [{"path": r[0], "title": r[1], "excerpt": r[2]} for r in rows]


# ─── AI Chat ──────────────────────────────────────────────────────────────────

class ChatRequest(BaseModel):
    message: str
    history: list = []


@app.get("/api/ai/status")
def ai_status():
    try:
        import sys as _sys
        _sys.path.insert(0, str(ROOT / "scripts"))
        from arcana_ai_agent import ArcanaAgent
        alive = ArcanaAgent.is_available()
        models = ArcanaAgent.models() if alive else []
        return {"available": alive, "models": models}
    except Exception as e:
        return {"available": False, "error": str(e)}


@app.post("/api/ai/chat")
async def ai_chat(req: ChatRequest):
    try:
        import sys as _sys
        _sys.path.insert(0, str(ROOT / "scripts"))
        from arcana_ai_agent import ArcanaAgent
        agent = ArcanaAgent()
        if not agent.is_available():
            return JSONResponse(
                status_code=503,
                content={"error": "Ollama is not running. Start it with: ollama serve"}
            )
        response, updated_history = agent.respond(req.history, req.message)
        return {"response": response, "history": updated_history}
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.backend.main:app", host="0.0.0.0", port=8000, reload=True)
