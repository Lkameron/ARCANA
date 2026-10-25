# SCRIPT: 06_lncrna_annotation.R
#
# Inputs:
#   - results/lncrna/tables/lncRNA_cross_comparison_summary.tsv
#   - results/lncrna/tables/lncRNA_points_of_interest_T2D_vs_ND_for_networks.tsv
#   - results/de/tables/DESeq2_T2D_vs_ND_lncRNA_only.tsv          (optional)
#   - results/de/tables/DESeq2_IGT_vs_ND_lncRNA_only.tsv          (optional)
#   - results/de/tables/DESeq2_T2D_vs_IGT_lncRNA_only.tsv         (optional)
#   - results/wgcna/modules/wgcna_gene_module_assignments.csv
#   - results/wgcna/module_trait/wgcna_module_status_associations_summary.csv
#   - results/wgcna/hubgenes/tables/wgcna_hub_genes_ranked_focus_modules.csv     (optional)
#   - results/wgcna/hubgenes/tables/wgcna_hub_lncRNAs_ranked_focus_modules.csv   (optional)
#   - data/metadata/gene_annotation.tsv                                  (optional)
#
# Outputs:
#   - results/lncrna/network/tables/lncRNA_WGCNA_module_integration.tsv
#   - results/lncrna/network/tables/lncRNA_module_hub_summary.tsv
#   - results/lncrna/network/tables/lncRNA_candidates_for_network_visualisation.tsv
#   - results/lncrna/network/plots/lncRNA_module_priority_heatmap.png
#   - results/lncrna/network/plots/lncRNA_module_vs_priority_scatter.png
#   - GSE164416_lncRNA_network_annotation_summary.xlsx

# setup

required_cran <- c(
  "data.table",
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "ggplot2",
  "ggrepel",
  "tidyr",
  "forcats",
  "writexl"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(data.table)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(forcats)
library(writexl)

options(stringsAsFactors = FALSE)

cat("========================================\n")
cat("06 - lncRNA annotation and WGCNA integration\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

# Define paths

project_dir <- getwd()

lnc_results_dir  <- file.path(project_dir, "results", "lncrna")
lnc_tables_dir   <- file.path(lnc_results_dir, "tables")
lnc_plots_dir    <- file.path(lnc_results_dir, "plots")

lnc_cross_summary_file <- file.path(lnc_tables_dir, "lncRNA_cross_comparison_summary.tsv")
lnc_points_file        <- file.path(lnc_tables_dir, "lncRNA_points_of_interest_T2D_vs_ND_for_networks.tsv")

de_tables_dir <- file.path(project_dir, "results", "de", "tables")
de_t2d_vs_nd_file  <- file.path(de_tables_dir, "DESeq2_T2D_vs_ND_lncRNA_only.tsv")
de_igt_vs_nd_file  <- file.path(de_tables_dir, "DESeq2_IGT_vs_ND_lncRNA_only.tsv")
de_t2d_vs_igt_file <- file.path(de_tables_dir, "DESeq2_T2D_vs_IGT_lncRNA_only.tsv")

wgcna_dir          <- file.path(project_dir, "results", "wgcna")
wgcna_tables_dir   <- file.path(wgcna_dir, "tables")
wgcna_modules_dir  <- file.path(wgcna_dir, "modules")
wgcna_hub_dir      <- file.path(wgcna_dir, "hubgenes", "tables")
wgcna_mt_dir       <- file.path(wgcna_dir, "module_trait")

module_assignments_file  <- file.path(wgcna_modules_dir, "wgcna_gene_module_assignments.csv")
module_status_assoc_file <- file.path(wgcna_mt_dir, "wgcna_module_status_associations_summary.csv")
hub_all_genes_file       <- file.path(wgcna_hub_dir, "wgcna_hub_genes_ranked_focus_modules.csv")
hub_lnc_file             <- file.path(wgcna_hub_dir, "wgcna_hub_lncRNAs_ranked_focus_modules.csv")

metadata_dir    <- file.path(project_dir, "data", "metadata")
gene_annot_file <- file.path(metadata_dir, "gene_annotation.tsv")

lnc_network_dir        <- file.path(lnc_results_dir, "network")
lnc_network_tables_dir <- file.path(lnc_network_dir, "tables")
lnc_network_plots_dir  <- file.path(lnc_network_dir, "plots")

dir.create(lnc_network_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(lnc_network_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lnc_network_plots_dir,  recursive = TRUE, showWarnings = FALSE)

# Plot settings

status_colors <- c(
  "ND"  = "#012345",
  "IGT" = "#345678",
  "T2D" = "#8b0000"
)

neutral_dark <- "#2B1B17"

theme_set(
  theme_bw(base_size = 12, base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
)

# Helper functions

check_file_exists_nonempty <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Required file not found: ", path)
    return(FALSE)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) {
    if (required) stop("File exists but is empty: ", path)
    return(FALSE)
  }
  TRUE
}

safe_read_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE)
}

pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

rescale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) return(rep(NA_real_, length(x)))
  if (diff(rng) == 0) return(ifelse(is.na(x), NA_real_, 1))
  (x - rng[1]) / (rng[2] - rng[1])
}

ensure_column <- function(df, colname, default = NA) {
  if (!colname %in% colnames(df)) {
    df[[colname]] <- default
  }
  df
}

# Check required inputs

cat("Checking required input files...\n")

required_files <- c(
  lnc_cross_summary_file,
  lnc_points_file,
  module_assignments_file,
  module_status_assoc_file
)

for (f in required_files) {
  check_file_exists_nonempty(f, required = TRUE)
}

has_gene_annot    <- check_file_exists_nonempty(gene_annot_file, required = FALSE)
has_de_t2d_vs_nd  <- check_file_exists_nonempty(de_t2d_vs_nd_file, required = FALSE)
has_de_igt_vs_nd  <- check_file_exists_nonempty(de_igt_vs_nd_file, required = FALSE)
has_de_t2d_vs_igt <- check_file_exists_nonempty(de_t2d_vs_igt_file, required = FALSE)
has_hub_all_genes <- check_file_exists_nonempty(hub_all_genes_file, required = FALSE)
has_hub_lnc       <- check_file_exists_nonempty(hub_lnc_file, required = FALSE)

cat("Input check complete.\n\n")

# Load core data

cat("Loading lncRNA prioritisation and WGCNA outputs...\n")

lnc_cross          <- safe_read_tsv(lnc_cross_summary_file)
lnc_points         <- safe_read_tsv(lnc_points_file)
module_assignments <- safe_read_csv(module_assignments_file)
module_status_assoc <- safe_read_csv(module_status_assoc_file)

gene_annot <- NULL
if (has_gene_annot) {
  gene_annot <- safe_read_tsv(gene_annot_file)
}

de_lnc_all <- list()
if (has_de_t2d_vs_nd)  de_lnc_all[["T2D_vs_ND"]]  <- safe_read_tsv(de_t2d_vs_nd_file)
if (has_de_igt_vs_nd)  de_lnc_all[["IGT_vs_ND"]]  <- safe_read_tsv(de_igt_vs_nd_file)
if (has_de_t2d_vs_igt) de_lnc_all[["T2D_vs_IGT"]] <- safe_read_tsv(de_t2d_vs_igt_file)

hub_all_genes <- if (has_hub_all_genes) safe_read_csv(hub_all_genes_file) else NULL
hub_lnc       <- if (has_hub_lnc)       safe_read_csv(hub_lnc_file)       else NULL

cat("lncRNA cross-comparison rows: ", nrow(lnc_cross), "\n", sep = "")
cat("lncRNA T2D vs ND points-of-interest rows: ", nrow(lnc_points), "\n", sep = "")
cat("Module assignments rows: ", nrow(module_assignments), "\n", sep = "")
cat("Module status association rows: ", nrow(module_status_assoc), "\n\n", sep = "")

# Harmonise identifiers

cat("Harmonising gene identifiers...\n")

safe_rename <- function(df, new_name, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0 || hit[1] == new_name) return(df)
  dplyr::rename(df, !!new_name := dplyr::all_of(hit[1]))
}

lnc_cross <- lnc_cross %>%
  safe_rename("gene_id", c("geneid", "gene_id")) %>%
  mutate(gene_id = as.character(gene_id))

lnc_points <- lnc_points %>%
  safe_rename("gene_id", c("geneid", "gene_id")) %>%
  mutate(gene_id = as.character(gene_id))

module_assignments <- module_assignments %>%
  safe_rename("gene_id",      c("gene_id", "geneid", "original_gene_id")) %>%
  safe_rename("module_color", c("module_color", "module"))

if (!"gene_id" %in% colnames(module_assignments)) {
  stop("Module assignments must contain a gene_id (or equivalent) column.")
}
if (!"module_color" %in% colnames(module_assignments)) {
  stop("Module assignments must contain a module_color (or equivalent) column.")
}

module_assignments <- module_assignments %>%
  mutate(
    gene_id = as.character(gene_id),
    module_color = as.character(module_color)
  )

if (!is.null(gene_annot)) {
  annot_gene_col <- pick_first_existing(gene_annot, c("original_gene_id", "gene_id", "geneid"))
  if (!is.na(annot_gene_col)) {
    gene_annot <- safe_rename(gene_annot, "gene_id", c(annot_gene_col))
    keep_cols <- intersect(c("gene_id", "gene_symbol", "gene_name", "biotype", "ensembl_gene_id"),
                           colnames(gene_annot))
    gene_annot <- gene_annot[, keep_cols, drop = FALSE]
    gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), ]
    gene_annot$gene_id <- as.character(gene_annot$gene_id)
  } else {
    warning("No usable gene ID column found in gene_annotation.tsv; skipping annotation join.")
    gene_annot <- NULL
  }
}

if (!is.null(hub_all_genes)) {
  hub_all_genes <- hub_all_genes %>%
    safe_rename("gene_id",      c("gene_id", "geneid")) %>%
    safe_rename("module_color", c("module_color", "module")) %>%
    mutate(
      gene_id = as.character(gene_id),
      module_color = as.character(module_color)
    )
}

if (!is.null(hub_lnc)) {
  hub_lnc <- hub_lnc %>%
    safe_rename("gene_id",      c("gene_id", "geneid")) %>%
    safe_rename("module_color", c("module_color", "module")) %>%
    mutate(
      gene_id = as.character(gene_id),
      module_color = as.character(module_color)
    )
}

# Summarise module–disease associations

cat("Summarising module–disease associations...\n")

required_assoc_cols <- c("module", "trait", "correlation", "p_value")
if (!all(required_assoc_cols %in% colnames(module_status_assoc))) {
  stop(
    "Module status association file must contain columns: ",
    paste(required_assoc_cols, collapse = ", ")
  )
}

module_status_long <- module_status_assoc %>%
  mutate(
    module_color = gsub("^ME", "", module),
    trait_clean = trait
  ) %>%
  filter(trait_clean %in% c("status_ND", "status_IGT", "status_T2D")) %>%
  group_by(module_color) %>%
  summarise(
    best_trait = trait_clean[which.min(p_value)],
    best_cor   = correlation[which.min(p_value)],
    best_p     = min(p_value, na.rm = TRUE),
    .groups    = "drop"
  )

# Build lncRNA–module integration table

cat("Building lncRNA–WGCNA module integration table...\n")

lnc_cross <- lnc_cross %>%
  mutate(
    biotype   = coalesce(biotype, ""),
    is_lncRNA = str_detect(
      tolower(biotype),
      "lncrna|linc|antisense|noncoding|non-coding|processed_transcript"
    )
  )

lnc_cross_lnc_only <- lnc_cross %>%
  filter(is_lncRNA)

lnc_points <- lnc_points %>%
  mutate(is_lncRNA = TRUE)

lnc_module_tbl <- lnc_cross_lnc_only %>%
  full_join(
    module_assignments %>%
      select(
        gene_id,
        module_color,
        any_of(c("gene_symbol", "genesymbol", "gene_name", "genename", "biotype"))
      ),
    by = "gene_id"
  )

if (!is.null(gene_annot)) {
  annot_cols <- c(
    "gene_id",
    intersect(colnames(gene_annot), c("genesymbol", "gene_symbol")),
    intersect(colnames(gene_annot), c("genename", "gene_name")),
    intersect(colnames(gene_annot), c("biotype", "gene_biotype", "gene_type"))
  )
  annot_cols <- unique(annot_cols[annot_cols %in% colnames(gene_annot)])
  gene_annot2 <- gene_annot[, annot_cols, drop = FALSE]

  lnc_module_tbl <- lnc_module_tbl %>%
    left_join(gene_annot2, by = "gene_id", suffix = c("_lnc", "_annot")) %>%
    mutate(
      gene_symbol = coalesce(
        if ("gene_symbol_lnc"   %in% colnames(.)) gene_symbol_lnc   else NA_character_,
        if ("genesymbol_lnc"    %in% colnames(.)) genesymbol_lnc    else NA_character_,
        if ("gene_symbol_annot" %in% colnames(.)) gene_symbol_annot else NA_character_,
        if ("genesymbol_annot"  %in% colnames(.)) genesymbol_annot  else NA_character_
      ),
      gene_name = coalesce(
        if ("gene_name_lnc"   %in% colnames(.)) gene_name_lnc   else NA_character_,
        if ("genename_lnc"    %in% colnames(.)) genename_lnc    else NA_character_,
        if ("gene_name_annot" %in% colnames(.)) gene_name_annot else NA_character_,
        if ("genename_annot"  %in% colnames(.)) genename_annot  else NA_character_
      ),
      biotype = coalesce(
        if ("biotype_lnc"   %in% colnames(.)) biotype_lnc   else NA_character_,
        if ("biotype_annot" %in% colnames(.)) biotype_annot else NA_character_,
        if ("biotype"       %in% colnames(.)) biotype       else NA_character_
      )
    )
}

lnc_module_tbl <- ensure_column(lnc_module_tbl, "gene_symbol", NA_character_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "gene_name", NA_character_)

lnc_module_tbl <- lnc_module_tbl %>%
  left_join(module_status_long, by = "module_color")

lnc_module_tbl <- lnc_module_tbl %>%
  left_join(
    lnc_points %>%
      select(
        gene_id,
        any_of(c("prioritisation_rank", "prioritisationrank")),
        any_of(c("network_priority_flag", "networkpriorityflag")),
        any_of(c("comparisonspresent", "comparisons_present")),
        any_of(c("recurrencescore")),
        any_of(c("minpadj", "min_padj")),
        any_of(c("prioritisationscore", "prioritisation_score")),
        any_of(c("bestrank"))
      ),
    by = "gene_id"
  )

colnames(lnc_module_tbl) <- gsub("prioritisationrank",  "prioritisation_rank",  colnames(lnc_module_tbl))
colnames(lnc_module_tbl) <- gsub("networkpriorityflag", "network_priority_flag", colnames(lnc_module_tbl))
colnames(lnc_module_tbl) <- gsub("min_padj",            "minpadj",               colnames(lnc_module_tbl))
colnames(lnc_module_tbl) <- gsub("comparisons_present", "comparisonspresent",    colnames(lnc_module_tbl))
colnames(lnc_module_tbl) <- gsub("prioritisation_score","prioritisationscore",   colnames(lnc_module_tbl))

# Successive joins against lnc_cross / module_status_long / lnc_points can each
# introduce their own copy of the same field, which dplyr then suffixes as
# ".x"/".y" instead of merging. Collapse every suffixed/alternate spelling of
# each scoring input back into a single canonical column so downstream
# rescale01() calls see the real joined values rather than an empty stand-in.
resolve_dup_col <- function(df, canonical, alt_names = character()) {
  candidates <- unique(c(canonical, paste0(canonical, c(".x", ".y")), alt_names))
  candidates <- candidates[candidates %in% colnames(df)]
  if (length(candidates) == 0) return(df)
  vals <- do.call(dplyr::coalesce, lapply(candidates, function(cn) df[[cn]]))
  df <- df[, !(colnames(df) %in% setdiff(candidates, canonical)), drop = FALSE]
  df[[canonical]] <- vals
  df
}

lnc_module_tbl <- lnc_module_tbl %>%
  resolve_dup_col("comparisonspresent",  c("comparisons_present")) %>%
  resolve_dup_col("recurrencescore",     c("recurrence_score")) %>%
  resolve_dup_col("bestrank",            c("best_rank")) %>%
  resolve_dup_col("prioritisationscore", c("prioritisation_score")) %>%
  resolve_dup_col("minpadj",             c("min_padj")) %>%
  resolve_dup_col("prioritisation_rank") %>%
  resolve_dup_col("network_priority_flag")

lnc_module_tbl <- ensure_column(lnc_module_tbl, "prioritisationscore", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "prioritisation_rank", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "bestrank", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "recurrencescore", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "comparisonspresent", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "minpadj", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "network_priority_flag", NA)

lnc_module_tbl <- lnc_module_tbl %>%
  mutate(
    prioritisationscore = suppressWarnings(as.numeric(prioritisationscore)),
    prioritisation_rank = suppressWarnings(as.numeric(prioritisation_rank)),
    bestrank            = suppressWarnings(as.numeric(bestrank)),
    recurrencescore     = suppressWarnings(as.numeric(recurrencescore)),
    comparisonspresent  = suppressWarnings(as.numeric(comparisonspresent)),
    minpadj             = suppressWarnings(as.numeric(minpadj))
  ) %>%
  mutate(
    prioritisation_rank = coalesce(prioritisation_rank, bestrank),
    recurrencescore     = coalesce(recurrencescore, comparisonspresent)
  )

# Merge hub metrics too (maybe)

cat("Merging hub metrics where available...\n")

if (!is.null(hub_all_genes)) {
  hub_subset <- hub_all_genes %>%
    select(
      gene_id,
      module_color,
      any_of(c(
        "own_kME", "abs_own_kME", "priority_gs",
        "hub_score_all_genes", "hub_score_lncRNA",
        "hub_rank_all_genes"
      ))
    )

  lnc_module_tbl <- lnc_module_tbl %>%
    left_join(hub_subset, by = c("gene_id", "module_color"))
}

if (!is.null(hub_lnc)) {
  hub_lnc_subset <- hub_lnc %>%
    select(
      gene_id,
      module_color,
      any_of(c("hub_rank_lncRNA"))
    )

  lnc_module_tbl <- lnc_module_tbl %>%
    left_join(hub_lnc_subset, by = c("gene_id", "module_color"))
}

lnc_module_tbl <- ensure_column(lnc_module_tbl, "own_kME", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "abs_own_kME", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "priority_gs", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "hub_score_lncRNA", NA_real_)
lnc_module_tbl <- ensure_column(lnc_module_tbl, "hub_rank_lncRNA", NA_real_)

# Score lncRNAs

cat("Scoring lncRNAs for network relevance...\n")

lnc_module_tbl <- lnc_module_tbl %>%
  mutate(
    score_lnc_priority = rescale01(prioritisationscore),
    score_recurrence   = rescale01(recurrencescore),
    score_module_cor   = rescale01(abs(best_cor)),
    score_hub_central  = rescale01(abs_own_kME),
    score_gs           = rescale01(priority_gs)
  ) %>%
  mutate(
    network_integration_score =
      0.30 * score_lnc_priority +
      0.20 * score_recurrence +
      0.20 * score_module_cor +
      0.15 * score_hub_central +
      0.15 * score_gs
  ) %>%
  arrange(
    desc(network_integration_score),
    desc(score_lnc_priority),
    desc(score_hub_central)
  )

# Candidate sets

cat("Deriving lncRNA candidate sets for network visualisation...\n")

lnc_candidates_network <- lnc_module_tbl %>%
  filter(
    !is.na(module_color),
    !is.na(network_integration_score)
  ) %>%
  group_by(module_color) %>%
  arrange(
    desc(network_integration_score),
    desc(score_lnc_priority),
    desc(score_hub_central),
    .by_group = TRUE
  ) %>%
  mutate(
    module_rank_network = row_number()
  ) %>%
  ungroup()

# Top lncRNAs per module (up to 10 per module) using row_number() instead of slice_head()
lnc_top_per_module <- lnc_candidates_network %>%
  group_by(module_color) %>%
  filter(row_number() <= 10L) %>%
  ungroup()

# Export tables
cat("Exporting lncRNA–WGCNA integration tables...\n")

readr::write_tsv(
  lnc_module_tbl,
  file.path(lnc_network_tables_dir, "lncRNA_WGCNA_module_integration.tsv")
)

readr::write_tsv(
  lnc_top_per_module,
  file.path(lnc_network_tables_dir, "lncRNA_module_hub_summary.tsv")
)

# coffee break ;)
network_export <- lnc_candidates_network %>%
  select(
    gene_id,
    gene_symbol,
    gene_name,
    biotype,
    module_color,
    best_trait,
    best_cor,
    best_p,
    prioritisationscore,
    prioritisation_rank,
    recurrencescore,
    comparisonspresent,
    network_priority_flag,
    own_kME,
    abs_own_kME,
    priority_gs,
    hub_score_lncRNA,
    hub_rank_lncRNA,
    network_integration_score
  )

readr::write_tsv(
  network_export,
  file.path(lnc_network_tables_dir, "lncRNA_candidates_for_network_visualisation.tsv")
)

excel_output <- list(
  lncRNA_module_integration = lnc_module_tbl,
  lncRNA_module_hub_summary = lnc_top_per_module,
  lncRNA_network_candidates = network_export
)

writexl::write_xlsx(
  excel_output,
  path = file.path(project_dir, "GSE164416_lncRNA_network_annotation_summary.xlsx")
)

# Plots

cat("Creating lncRNA–module integration plots...\n")

heatmap_df <- lnc_top_per_module %>%
  mutate(
    gene_label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id)
  )

if (nrow(heatmap_df) > 0) {
  p_heat <- ggplot(
    heatmap_df,
    aes(x = module_color, y = fct_reorder(gene_label, network_integration_score))
  ) +
    geom_tile(aes(fill = network_integration_score), color = "white") +
    scale_fill_gradient(
      low = "#E6F0FA",
      high = "#8b0000",
      name = "Integration\nscore"
    ) +
    labs(
      title = "Top lncRNA candidates per WGCNA module",
      subtitle = "Scores integrate lncRNA prioritisation with module and hub context",
      x = "Module colour",
      y = "lncRNA"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = file.path(lnc_network_plots_dir, "lncRNA_module_priority_heatmap.png"),
    plot = p_heat,
    width = 9.5,
    height = 8,
    dpi = 300
  )
}

scatter_df <- lnc_candidates_network %>%
  mutate(
    label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id)
  )

if (nrow(scatter_df) > 0) {
  p_scatter <- ggplot(
    scatter_df,
    aes(x = best_cor, y = network_integration_score, color = best_trait)
  ) +
    geom_point(alpha = 0.75, size = 2.1) +
    geom_smooth(method = "lm", se = FALSE, color = neutral_dark, linewidth = 0.8) +
    scale_color_manual(
      values = c(
        "status_ND"  = status_colors["ND"],
        "status_IGT" = status_colors["IGT"],
        "status_T2D" = status_colors["T2D"]
      ),
      na.translate = FALSE
    ) +
    labs(
      title = "lncRNA network integration vs module–trait association",
      subtitle = "Higher scores in strongly T2D-associated modules suggest T2D-relevant lncRNA hubs",
      x = "Module–trait correlation (best disease-status association)",
      y = "lncRNA network integration score",
      color = "Trait"
    )

  ggsave(
    filename = file.path(lnc_network_plots_dir, "lncRNA_module_vs_priority_scatter.png"),
    plot = p_scatter,
    width = 9.5,
    height = 7.5,
    dpi = 300
  )
}

# SECTION 15: Final summary

cat("\n========================================\n")
cat("06 COMPLETE - lncRNA annotation and network integration\n")
cat("========================================\n")
cat("Total rows in lncRNA cross-comparison table: ", nrow(lnc_cross), "\n", sep = "")
cat("Total lncRNAs with module annotation: ", sum(!is.na(lnc_module_tbl$module_color)), "\n", sep = "")
cat("lncRNAs selected for network export: ", nrow(network_export), "\n", sep = "")
cat("Outputs written to:\n")
cat(" - ", file.path(lnc_network_tables_dir, "lncRNA_WGCNA_module_integration.tsv"), "\n", sep = "")
cat(" - ", file.path(lnc_network_tables_dir, "lncRNA_module_hub_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(lnc_network_tables_dir, "lncRNA_candidates_for_network_visualisation.tsv"), "\n", sep = "")
cat(" - ", file.path(lnc_network_plots_dir, "lncRNA_module_priority_heatmap.png"), "\n", sep = "")
cat(" - ", file.path(lnc_network_plots_dir, "lncRNA_module_vs_priority_scatter.png"), "\n", sep = "")
cat(" - ", file.path(project_dir, "GSE164416_lncRNA_network_annotation_summary.xlsx"), "\n", sep = "")
cat("========================================\n")
cat("Next step: run 07_enrichment_analysis.R on key modules and lncRNA-enriched hubs.\n")
cat("========================================\n")