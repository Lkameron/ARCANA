# config_utils.R ─ Source this at the top of every ARCANA R script.
# Reads project settings from environment variables set by the backend
# (which in turn reads them from config/project.yaml).
#
# Usage:
#   source(file.path(Sys.getenv("ARCANA_SCRIPTS_DIR", "scripts"), "config_utils.R"))
#
# All settings live in the list: arcana_cfg

suppressPackageStartupMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", repos = "https://cloud.r-project.org")
  library(jsonlite)
})

.parse_list <- function(x) {
  v <- trimws(unlist(strsplit(x, ",")))
  v[nchar(v) > 0]
}

.nullable_int <- function(x) {
  if (is.null(x) || nchar(trimws(x)) == 0) NULL else as.integer(x)
}

arcana_cfg <- list(
  # ── Project ──────────────────────────────────────────────────────────────────
  project_name      = Sys.getenv("ARCANA_PROJECT_NAME", "ARCANA Project"),
  dataset_id        = Sys.getenv("ARCANA_DATASET_ID", ""),
  description       = Sys.getenv("ARCANA_DESCRIPTION", ""),

  # ── Species ───────────────────────────────────────────────────────────────────
  species           = Sys.getenv("ARCANA_SPECIES", "human"),
  genome            = Sys.getenv("ARCANA_GENOME", "GRCh38"),
  ensembl_version   = as.integer(Sys.getenv("ARCANA_ENSEMBL_VERSION", "109")),

  # ── Input paths ───────────────────────────────────────────────────────────────
  counts_matrix     = Sys.getenv("ARCANA_COUNTS_MATRIX"),
  metadata_path     = Sys.getenv("ARCANA_METADATA"),
  results_dir       = Sys.getenv("ARCANA_RESULTS_DIR", "results"),

  # ── Metadata columns ──────────────────────────────────────────────────────────
  sample_id_col     = Sys.getenv("ARCANA_SAMPLE_ID_COL", "sample_id"),
  condition_col     = Sys.getenv("ARCANA_CONDITION_COL", "condition"),

  # ── Comparisons (list of named lists with name/group_a/group_b/label) ────────
  comparisons       = {
    raw <- Sys.getenv("ARCANA_COMPARISONS_JSON", "[]")
    tryCatch(jsonlite::fromJSON(raw, simplifyDataFrame = FALSE), error = function(e) list())
  },
  comparison_names  = .parse_list(Sys.getenv("ARCANA_COMPARISONS", "")),

  # ── Traits ───────────────────────────────────────────────────────────────────
  traits_continuous = .parse_list(Sys.getenv("ARCANA_TRAITS_CONTINUOUS", "")),
  traits_binary     = .parse_list(Sys.getenv("ARCANA_TRAITS_BINARY", "")),

  # ── DE ───────────────────────────────────────────────────────────────────────
  padj              = as.numeric(Sys.getenv("ARCANA_PADJ", "0.05")),
  log2fc            = as.numeric(Sys.getenv("ARCANA_LOG2FC", "1.0")),
  de_covariates     = .parse_list(Sys.getenv("ARCANA_DE_COVARIATES", "")),

  # ── WGCNA ────────────────────────────────────────────────────────────────────
  soft_power        = .nullable_int(Sys.getenv("ARCANA_WGCNA_SOFT_POWER", "")),
  min_module_size   = as.integer(Sys.getenv("ARCANA_WGCNA_MIN_MODULE", "30")),
  merge_cut_height  = as.numeric(Sys.getenv("ARCANA_WGCNA_MERGE_CUT", "0.25")),
  n_focus_modules   = as.integer(Sys.getenv("ARCANA_WGCNA_N_FOCUS", "5")),

  # ── lncRNA ───────────────────────────────────────────────────────────────────
  lncrna_enabled    = Sys.getenv("ARCANA_LNCRNA_ENABLED", "TRUE") == "TRUE",
  lncrna_biotypes   = .parse_list(Sys.getenv("ARCANA_LNCRNA_BIOTYPES", "lncRNA")),

  # ── Figures ──────────────────────────────────────────────────────────────────
  figures_dpi       = as.integer(Sys.getenv("ARCANA_FIGURES_DPI", "300"))
)

# ── Derived paths ─────────────────────────────────────────────────────────────
arcana_cfg$de_dir        <- file.path(arcana_cfg$results_dir, "de")
arcana_cfg$wgcna_dir     <- file.path(arcana_cfg$results_dir, "wgcna")
arcana_cfg$ml_dir        <- file.path(arcana_cfg$results_dir, "ml")
arcana_cfg$enrichment_dir <- file.path(arcana_cfg$results_dir, "enrichment")
arcana_cfg$figures_dir   <- file.path(arcana_cfg$results_dir, "figures")
arcana_cfg$reports_dir   <- file.path(arcana_cfg$results_dir, "reports")

# ── Species-to-OrgDb mapping (for enrichment) ─────────────────────────────────
arcana_cfg$orgdb <- switch(arcana_cfg$species,
  "human" = "org.Hs.eg.db",
  "mouse" = "org.Mm.eg.db",
  "rat"   = "org.Rn.eg.db",
  "org.Hs.eg.db"   # default to human
)

# ── Create output directories ─────────────────────────────────────────────────
for (.d in c(arcana_cfg$de_dir, arcana_cfg$wgcna_dir, arcana_cfg$ml_dir,
             arcana_cfg$enrichment_dir, arcana_cfg$figures_dir, arcana_cfg$reports_dir)) {
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}

message("[ARCANA] Config loaded — project: ", arcana_cfg$project_name,
        " | comparisons: ", paste(arcana_cfg$comparison_names, collapse = ", "))
