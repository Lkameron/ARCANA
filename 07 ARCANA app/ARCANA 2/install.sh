#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# ARCANA — install.sh
# Installs all dependencies for running ARCANA locally (without Docker).
# Run once before starting the platform for the first time.
#
# Usage: bash install.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

PYTHON=${PYTHON:-python3.11}
ARCANA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""
echo "╔══════════════════════════════════════╗"
echo "║        ARCANA — Installer            ║"
echo "╚══════════════════════════════════════╝"
echo ""

# ── 1. Check Python 3.11 ──────────────────────────────────────────────────────
echo "[1/5] Checking Python 3.11..."
if ! command -v "$PYTHON" &>/dev/null; then
  echo "  ERROR: $PYTHON not found."
  echo "  Install it from https://www.python.org/downloads/ or via pyenv."
  exit 1
fi
PY_VER=$("$PYTHON" --version 2>&1)
echo "  Found: $PY_VER"

# ── 2. Create virtual environment ─────────────────────────────────────────────
echo "[2/5] Creating virtual environment..."
if [ ! -d "$ARCANA_DIR/venv" ]; then
  "$PYTHON" -m venv "$ARCANA_DIR/venv"
  echo "  Created: $ARCANA_DIR/venv"
else
  echo "  Exists: $ARCANA_DIR/venv (skipping)"
fi

source "$ARCANA_DIR/venv/bin/activate"
pip install --quiet --upgrade pip

# ── 3. Install Python dependencies ────────────────────────────────────────────
echo "[3/5] Installing Python packages..."
pip install --quiet -r "$ARCANA_DIR/requirements.txt"
echo "  Done."

# ── 4. Check R ───────────────────────────────────────────────────────────────
echo "[4/5] Checking R..."
if ! command -v Rscript &>/dev/null; then
  echo "  WARNING: R not found. Install from https://www.r-project.org/"
  echo "  R is required to run the analysis pipeline."
else
  R_VER=$(Rscript --version 2>&1 | head -1)
  echo "  Found: $R_VER"
  echo "  Installing R packages (this may take several minutes on first run)..."
  Rscript - <<'RSCRIPT'
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    pkgs_cran <- c("BiocManager","jsonlite","data.table","dplyr","readr",
                   "tibble","ggplot2","ggrepel","writexl","stringr","tidyr",
                   "forcats","yaml")
    to_install <- pkgs_cran[!vapply(pkgs_cran, requireNamespace, logical(1), quietly = TRUE)]
    if (length(to_install)) install.packages(to_install)

    pkgs_bioc <- c("DESeq2","WGCNA","clusterProfiler","EnhancedVolcano",
                   "ComplexHeatmap","circlize","org.Hs.eg.db","org.Mm.eg.db",
                   "enrichplot","AnnotationDbi")
    to_install_bioc <- pkgs_bioc[!vapply(pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]
    if (length(to_install_bioc)) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
    cat("R packages: OK\n")
RSCRIPT
fi

# ── 5. Check / install Ollama ─────────────────────────────────────────────────
echo "[5/5] Checking Ollama (AI assistant)..."
if command -v ollama &>/dev/null; then
  echo "  Ollama is installed."
  echo "  Pulling llama3.2 model (if not already downloaded)..."
  ollama pull llama3.2 || echo "  (Start 'ollama serve' first, then re-run this step)"
else
  echo "  Ollama not found."
  if [[ "$(uname)" == "Darwin" ]]; then
    echo "  Download the Mac app from: https://ollama.com/download"
    echo "  Then run: ollama pull llama3.2"
  else
    echo "  Install with: curl -fsSL https://ollama.com/install.sh | sh"
    echo "  Then run: ollama pull llama3.2"
  fi
fi

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════"
echo " Installation complete."
echo ""
echo " To start ARCANA:"
echo "   source venv/bin/activate"
echo "   uvicorn app.backend.main:app --port 8000"
echo "   Then open: http://localhost:8000"
echo ""
echo " Edit config/project.yaml before your first analysis run."
echo "═══════════════════════════════════════════════════════"
echo ""
