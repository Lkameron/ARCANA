# 05C_hub_gene_prioritisation.R
#
# Project:
#   Pancreatic-Islet-T2D-Coexpression
# SECTION 1: Package setup

required_cran <- c(
  "data.table",
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "ggplot2",
  "tidyr",
  "forcats"
)

required_bioc <- c("WGCNA")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(ggplot2)
library(tidyr)
library(forcats)
library(WGCNA)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()
allowWGCNAThreads()

cat("========================================\n")
cat("05C - Hub gene prioritisation\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

# Define paths

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else "run1"
cat("RUN_ID =", RUN_ID, "\n\n")

project_dir <- getwd()

# WGCNA inputs (TRAIN-only, from N1/N2 — nested-validation leakage fix)
input_expr_file <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables", "wgcna_expression_input_TRAIN.tsv")
input_traits_file <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables", "wgcna_traits_input_TRAIN.csv")
input_module_assignments <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "modules", "wgcna_gene_module_assignments.csv")
input_module_sizes <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "modules", "wgcna_module_sizes.csv")
input_05b_rdata <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "rdata", "05B_wgcna_network_construction.RData")

#  module-trait summary
input_module_trait_summary <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "module_trait", "wgcna_module_status_associations_summary.csv")

# external evidence (TRAIN-only DE fit, from N0 — closes the secondary leak
# where DE p-values/effect sizes computed on all 98 samples fed hub scoring)
input_de_t2d_vs_nd <- file.path(project_dir, "results", "nested_validation", RUN_ID, "de_train", "tables", "DESeq2_T2D_vs_ND_all_genes.tsv")
input_de_igt_vs_nd <- file.path(project_dir, "results", "nested_validation", RUN_ID, "de_train", "tables", "DESeq2_IGT_vs_ND_all_genes.tsv")
input_de_t2d_vs_igt <- file.path(project_dir, "results", "nested_validation", RUN_ID, "de_train", "tables", "DESeq2_T2D_vs_IGT_all_genes.tsv")
input_lnc_priority <- file.path(project_dir, "results", "lncrna", "tables", "lncRNA_points_of_interest_T2D_vs_ND_for_networks.tsv")

# Outputs
output_dir <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "hubgenes")
output_tables_dir <- file.path(output_dir, "tables")
output_plots_dir <- file.path(output_dir, "plots")
output_rdata_dir <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train", "rdata")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_rdata_dir, recursive = TRUE, showWarnings = FALSE)

# Plot styling

# Disease-status palette inherited from the rest of the project
status_colors <- c(
  "ND"  = "#012345",
  "IGT" = "#345678",
  "T2D" = "#8B0000"
)

neutral_dark <- "#2B1B17"

base_theme <- theme_bw(base_size = 12, base_family = "serif") +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88")
  )

# Helper functions

check_file <- function(path, required = TRUE) {
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

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE)
}

safe_read_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

safe_read_table <- function(path) {
  data.table::fread(path, data.table = FALSE)
}

pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

minmax_rescale <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(rep(NA_real_, length(x)))
  if (isTRUE(all.equal(rng[1], rng[2]))) return(rep(1, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

safe_cor <- function(x, y) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "pearson"))
}

safe_cor_test <- function(x, y) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 3) return(c(estimate = NA_real_, p.value = NA_real_))
  out <- tryCatch(cor.test(x[keep], y[keep], method = "pearson"), error = function(e) NULL)
  if (is.null(out)) return(c(estimate = NA_real_, p.value = NA_real_))
  c(estimate = unname(out$estimate), p.value = out$p.value)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

is_lnc_biotype <- function(x) {
  x <- tolower(as.character(x))
  stringr::str_detect(
    x,
    "lncrna|linc|antisense|sense_intronic|sense_overlapping|processed_transcript|non_coding|macro_lncRNA"
  )
}

# Check and load input files

required_files <- c(
  input_expr_file,
  input_traits_file,
  input_module_assignments,
  input_module_sizes,
  input_05b_rdata
)

for (f in required_files) {
  check_file(f, required = TRUE)
}

has_module_trait_summary <- check_file(input_module_trait_summary, required = FALSE)
has_de_t2d_vs_nd <- check_file(input_de_t2d_vs_nd, required = FALSE)
has_de_igt_vs_nd <- check_file(input_de_igt_vs_nd, required = FALSE)
has_de_t2d_vs_igt <- check_file(input_de_t2d_vs_igt, required = FALSE)
has_lnc_priority <- check_file(input_lnc_priority, required = FALSE)

cat("Input checks complete.\n\n")

cat("Loading WGCNA expression, trait, module-assignment, and module-size files...\n")

expr_df <- safe_read_table(input_expr_file)
traits_df <- safe_read_csv(input_traits_file)
module_assignments <- safe_read_csv(input_module_assignments)
module_sizes <- safe_read_csv(input_module_sizes)

load(input_05b_rdata)

cat("Expression rows: ", nrow(expr_df), " | columns: ", ncol(expr_df), "\n", sep = "")
cat("Trait rows: ", nrow(traits_df), " | columns: ", ncol(traits_df), "\n", sep = "")
cat("Module assignments rows: ", nrow(module_assignments), "\n", sep = "")
cat("Module sizes rows: ", nrow(module_sizes), "\n\n", sep = "")

# Prepare expression matrix and traits

if (!("sampleid" %in% colnames(expr_df))) {
  stop("The expression input file must contain a 'sampleid' column.")
}
if (!("sampleid" %in% colnames(traits_df))) {
  stop("The trait input file must contain a 'sampleid' column.")
}

datExpr <- expr_df %>%
  tibble::column_to_rownames("sampleid") %>%
  as.data.frame()

datExpr[] <- lapply(datExpr, as.numeric)
rownames(datExpr) <- as.character(rownames(datExpr))
traits_df$sampleid <- as.character(traits_df$sampleid)

common_samples <- intersect(rownames(datExpr), traits_df$sampleid)

if (length(common_samples) < 10) {
  stop("Too few overlapping samples between expression and trait files.")
}

datExpr <- datExpr[common_samples, , drop = FALSE]
traits_df <- traits_df %>%
  filter(sampleid %in% common_samples) %>%
  dplyr::slice(match(common_samples, sampleid))

stopifnot(all(rownames(datExpr) == traits_df$sampleid))

# Standardise core clinical variables
traits_df <- traits_df %>%
  mutate(
    status = factor(as.character(status), levels = c("ND", "IGT", "T2D")),
    status_ND  = ifelse(status == "ND", 1, 0),
    status_IGT = ifelse(status == "IGT", 1, 0),
    status_T2D = ifelse(status == "T2D", 1, 0),
    age = suppressWarnings(as.numeric(age)),
    bmi = suppressWarnings(as.numeric(bmi)),
    hba1c = suppressWarnings(as.numeric(hba1c)),
    glycaemia_2h = suppressWarnings(as.numeric(glycaemia_2h)),
    size_factor = suppressWarnings(as.numeric(size_factor))
  )

cat("Prepared WGCNA matrix: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes\n\n", sep = "")

# Choose modules of interest

default_focus_modules <- c("brown", "yellow", "green", "red", "black")

cat("Selecting modules of interest...\n")

focus_modules <- default_focus_modules
module_focus_reason <- "Default biologically interpretable mid-sized modules"

module_trait_summary <- NULL

if (has_module_trait_summary) {
  module_trait_summary <- safe_read_csv(input_module_trait_summary)

  required_cols <- c("module", "trait", "correlation", "p_value")
  if (all(required_cols %in% colnames(module_trait_summary))) {

    module_interest <- module_trait_summary %>%
      filter(trait %in% c("status_ND", "status_IGT", "status_T2D")) %>%
      mutate(module_clean = gsub("^ME", "", module)) %>%
      filter(module_clean != "grey") %>%
      arrange(p_value, desc(abs(correlation)))

    module_interest_unique <- module_interest %>%
      distinct(module_clean, .keep_all = TRUE)

    # Prefer modules that are both trait-associated and in the sensible
    overlap_modules <- intersect(default_focus_modules, module_interest_unique$module_clean)

    if (length(overlap_modules) >= 3) {
      focus_modules <- overlap_modules
      module_focus_reason <- "Overlap between trait-associated modules and interpretable mid-sized modules"
    } else if (nrow(module_interest_unique) > 0) {
      n_keep <- min(5L, nrow(module_interest_unique))
      focus_modules <- module_interest_unique$module_clean[seq_len(n_keep)]
      module_focus_reason <- "Top disease-status-associated modules from module-trait analysis"
    }
  }
}

focus_modules <- unique(focus_modules)
focus_modules <- focus_modules[focus_modules %in% module_sizes$module_color]

# Need to look at what lncRNAs are in the top modules and whether they have prioritisation evidence
# add in  the top lncRNA modules if they are not already included in the focus_modules list
lncrna_priority_file <- file.path(project_dir, "results", "lncrna", "tables",
                                   "lncRNA_cross_comparison_summary.tsv")
if (file.exists(lncrna_priority_file) && file.info(lncrna_priority_file)$size > 0) {
  lnc_tbl <- tryCatch(
    data.table::fread(lncrna_priority_file, data.table = FALSE),
    error = function(e) NULL
  )
  # Nested-validation fix: use the TRAIN-only module assignments (input_module_assignments),
  # not the original full-98-sample assignments, to stay fully leakage-free.
  wgcna_assign_file <- input_module_assignments
  if (!is.null(lnc_tbl) && file.exists(wgcna_assign_file)) {
    wgcna_assign <- tryCatch(
      data.table::fread(wgcna_assign_file, data.table = FALSE),
      error = function(e) NULL
    )
    if (!is.null(wgcna_assign)) {
      # Get top 20 lncRNA candidates by priority score
      score_col <- intersect(c("global_priority_score", "prioritisation_score", "mean_rank"),
                             colnames(lnc_tbl))[1]
      id_col    <- intersect(c("gene_id", "ensembl_gene_id"), colnames(lnc_tbl))[1]
      if (!is.na(score_col) && !is.na(id_col)) {
        if (score_col == "mean_rank") {
          top_lnc_ids <- lnc_tbl[order(lnc_tbl[[score_col]]), id_col][seq_len(min(20, nrow(lnc_tbl)))]
        } else {
          top_lnc_ids <- lnc_tbl[order(-lnc_tbl[[score_col]]), id_col][seq_len(min(20, nrow(lnc_tbl)))]
        }
        assign_id_col <- intersect(c("gene_id", "ensembl_gene_id"), colnames(wgcna_assign))[1]
        if (!is.na(assign_id_col)) {
          lnc_modules <- wgcna_assign[wgcna_assign[[assign_id_col]] %in% top_lnc_ids, "module_color"]
          lnc_modules <- unique(as.character(lnc_modules))
          lnc_modules <- lnc_modules[!is.na(lnc_modules) & lnc_modules != "grey" &
                                       lnc_modules %in% module_sizes$module_color]
          if (length(lnc_modules) > 0) {
            added <- setdiff(lnc_modules, focus_modules)
            if (length(added) > 0) {
              cat("Adding modules containing top lncRNA candidates: ",
                  paste(added, collapse = ", "), "\n", sep = "")
            }
            focus_modules <- unique(c(focus_modules, lnc_modules))
          }
        }
      }
    }
  }
}

if (length(focus_modules) == 0) {
  stop("No valid focus modules were identified.")
}

focus_module_table <- module_sizes %>%
  filter(module_color %in% focus_modules) %>%
  arrange(desc(n_genes)) %>%
  mutate(selection_reason = module_focus_reason)

readr::write_csv(
  focus_module_table,
  file.path(output_tables_dir, "focus_modules_selected_for_hub_analysis.csv")
)

cat("Focus modules selected: ", paste(focus_modules, collapse = ", "), "\n", sep = "")
cat("Selection rationale: ", module_focus_reason, "\n\n", sep = "")

# Prepare module assignment table

required_module_cols <- c("gene_id", "module_color")
missing_module_cols <- setdiff(required_module_cols, colnames(module_assignments))
if (length(missing_module_cols) > 0) {
  stop("Missing required columns in module assignment file: ",
       paste(missing_module_cols, collapse = ", "))
}

hub_tbl <- module_assignments %>%
  mutate(
    gene_id = as.character(gene_id),
    module_color = as.character(module_color)
  )

# Restrict to selected modules for main downstream prioritisation
hub_tbl_focus <- hub_tbl %>%
  filter(module_color %in% focus_modules)

cat("Genes in selected modules: ", nrow(hub_tbl_focus), "\n\n", sep = "")

# Extract module eigengenes from 05B RData

# The upstream 05B script may have stored module eigengenes under slightly
# different object names, so search for candidates.
candidate_me_names <- c("module_eigengenes", "MEs", "MEs_ordered", "MEs0")
me_object_name <- candidate_me_names[candidate_me_names %in% ls()][1]

if (is.na(me_object_name)) {
  stop("Could not find a module eigengene object in 05B_wgcna_network_construction.RData")
}

ME_mat <- get(me_object_name)
ME_mat <- as.data.frame(ME_mat)

if (nrow(ME_mat) != nrow(datExpr)) {
  # Try to align rows by sample names if I can 
  if (!is.null(rownames(ME_mat)) && all(rownames(datExpr) %in% rownames(ME_mat))) {
    ME_mat <- ME_mat[rownames(datExpr), , drop = FALSE]
  }
}

if (nrow(ME_mat) != nrow(datExpr)) {
  stop("Module eigengene matrix could not be aligned to the expression matrix.")
}

cat("Module eigengene object used: ", me_object_name, "\n", sep = "")
cat("Module eigengene dimensions: ", nrow(ME_mat), " samples x ", ncol(ME_mat), " eigengenes\n\n", sep = "")

# Calculate module membership (kME)

# Module membership (kME) = correlation of each gene with each module eigengene.
cat("Calculating module membership (kME) values...\n")

kME_mat <- as.data.frame(cor(datExpr, ME_mat, use = "pairwise.complete.obs", method = "pearson"))
colnames(kME_mat) <- paste0("kME_", colnames(kME_mat))
kME_mat$gene_id <- rownames(kME_mat)

# Calculate gene significance (GS)

# Gene significance is calculated here as correlation between each gene and
# disease status / clinical traits.
cat("Calculating gene significance (GS) values...\n")

trait_vars <- c("status_ND", "status_IGT", "status_T2D", "hba1c", "bmi", "glycaemia_2h")
trait_vars <- trait_vars[trait_vars %in% colnames(traits_df)]

gs_list <- list()

for (tr in trait_vars) {
  gs_cor <- apply(datExpr, 2, function(g) safe_cor(g, traits_df[[tr]]))
  gs_p <- apply(datExpr, 2, function(g) safe_cor_test(g, traits_df[[tr]])["p.value"])

  gs_tbl <- tibble(
    gene_id = names(gs_cor),
    !!paste0("GS_", tr) := as.numeric(gs_cor),
    !!paste0("GSabs_", tr) := abs(as.numeric(gs_cor)),
    !!paste0("GSp_", tr) := as.numeric(gs_p)
  )

  gs_list[[tr]] <- gs_tbl
}

gs_merged <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), gs_list)

# Merge kME and GS into focus table

hub_tbl_focus <- hub_tbl_focus %>%
  left_join(kME_mat, by = "gene_id") %>%
  left_join(gs_merged, by = "gene_id")

# SECTION 13: Determine each gene's own-module kME

# For each gene, extract the kME value corresponding to its own module.
get_own_kme <- function(module_color_value, row_values, kme_columns) {
  possible_names <- c(
    paste0("kME_ME", module_color_value),
    paste0("kME_", module_color_value),
    paste0("kME_ME", gsub("^ME", "", module_color_value))
  )
  hit <- possible_names[possible_names %in% kme_columns][1]
  if (is.na(hit)) return(NA_real_)
  as.numeric(row_values[[hit]])
}

kme_cols <- grep("^kME_", colnames(hub_tbl_focus), value = TRUE)

hub_tbl_focus$own_kME <- apply(
  hub_tbl_focus,
  1,
  function(r) get_own_kme(r[["module_color"]], as.list(r), kme_cols)
)

hub_tbl_focus <- hub_tbl_focus %>%
  mutate(abs_own_kME = abs(as.numeric(own_kME)))

# Add optional differential-expression evidence

merge_de_table <- function(main_tbl, de_path, prefix) {
  if (!file.exists(de_path)) return(main_tbl)

  de_tbl <- safe_read_tsv(de_path)
  gene_col <- pick_first_existing(de_tbl, c("gene_id", "GeneID", "gene", "id"))

  if (is.na(gene_col)) return(main_tbl)

  de_tbl <- de_tbl %>%
    rename(gene_id = all_of(gene_col))

  keep_cols <- c(
    "gene_id",
    intersect(c("log2FoldChange", "padj", "pvalue", "baseMean", "stat"), colnames(de_tbl))
  )

  de_tbl <- de_tbl[, keep_cols, drop = FALSE]
  non_gene_cols <- setdiff(colnames(de_tbl), "gene_id")
  colnames(de_tbl)[colnames(de_tbl) %in% non_gene_cols] <- paste0(prefix, "_", non_gene_cols)

  left_join(main_tbl, de_tbl, by = "gene_id")
}

cat("Adding differential-expression evidence where available...\n")

hub_tbl_focus <- merge_de_table(hub_tbl_focus, input_de_t2d_vs_nd, "T2D_vs_ND")
hub_tbl_focus <- merge_de_table(hub_tbl_focus, input_de_igt_vs_nd, "IGT_vs_ND")
hub_tbl_focus <- merge_de_table(hub_tbl_focus, input_de_t2d_vs_igt, "T2D_vs_IGT")

# Add optional lncRNA-prioritisation evidence

if (has_lnc_priority) {
  cat("Adding lncRNA prioritisation evidence...\n")

  lnc_priority_tbl <- safe_read_tsv(input_lnc_priority)
  lnc_gene_col <- pick_first_existing(lnc_priority_tbl, c("gene_id", "GeneID", "gene"))

  if (!is.na(lnc_gene_col)) {
    lnc_priority_tbl <- lnc_priority_tbl %>%
      rename(gene_id = all_of(lnc_gene_col))

    hub_tbl_focus <- left_join(hub_tbl_focus, lnc_priority_tbl, by = "gene_id")
  }
}

# Define composite hub-prioritisation scores

# The plan [12/05/2026]
# - For all genes: combine module membership, gene significance, DE significance,
#   and DE effect size.
# - For lncRNAs: also incorporate external lncRNA prioritisation score when
#   available.
cat("Constructing composite hub scores...\n")

# Prefer T2D association as the main disease-oriented GS measure if available.
if ("GSabs_status_T2D" %in% colnames(hub_tbl_focus)) {
  hub_tbl_focus$priority_gs <- hub_tbl_focus$GSabs_status_T2D
} else {
  gs_priority_cols <- intersect(
    c("GSabs_status_IGT", "GSabs_status_ND"),
    colnames(hub_tbl_focus)
  )
  if (length(gs_priority_cols) > 0) {
    hub_tbl_focus$priority_gs <- apply(
      hub_tbl_focus[, gs_priority_cols, drop = FALSE],
      1,
      function(x) {
        x <- as.numeric(x)
        x <- x[is.finite(x)]
        if (length(x) == 0) return(NA_real_)
        max(x)
      }
    )
  } else {
    hub_tbl_focus$priority_gs <- NA_real_
  }
}

if ("T2D_vs_ND_padj" %in% colnames(hub_tbl_focus)) {
  hub_tbl_focus$de_strength <- -log10(pmax(hub_tbl_focus$T2D_vs_ND_padj, 1e-300))
} else {
  hub_tbl_focus$de_strength <- NA_real_
}

if ("T2D_vs_ND_log2FoldChange" %in% colnames(hub_tbl_focus)) {
  hub_tbl_focus$de_lfc_abs <- abs(hub_tbl_focus$T2D_vs_ND_log2FoldChange)
} else {
  hub_tbl_focus$de_lfc_abs <- NA_real_
}

if ("prioritisation_score" %in% colnames(hub_tbl_focus)) {
  hub_tbl_focus$external_lnc_score <- suppressWarnings(as.numeric(hub_tbl_focus$prioritisation_score))
} else {
  hub_tbl_focus$external_lnc_score <- NA_real_
}

hub_tbl_focus <- hub_tbl_focus %>%
  mutate(
    score_kME = minmax_rescale(abs_own_kME),
    score_GS = minmax_rescale(priority_gs),
    score_DE = minmax_rescale(de_strength),
    score_LFC = minmax_rescale(de_lfc_abs),
    score_lnc_external = minmax_rescale(external_lnc_score)
  ) %>%
  mutate(
    hub_score_all_genes = rowMeans(
      cbind(score_kME, score_GS, score_DE, score_LFC),
      na.rm = TRUE
    ),
    hub_score_lncRNA = rowMeans(
      cbind(score_kME, score_GS, score_DE, score_LFC, score_lnc_external),
      na.rm = TRUE
    )
  )

hub_tbl_focus$hub_score_all_genes[!is.finite(hub_tbl_focus$hub_score_all_genes)] <- NA_real_
hub_tbl_focus$hub_score_lncRNA[!is.finite(hub_tbl_focus$hub_score_lncRNA)] <- NA_real_

# Rank hub genes within each selected module

cat("Ranking hub genes within each selected module...\n")

hub_ranked_all <- hub_tbl_focus %>%
  group_by(module_color) %>%
  arrange(desc(hub_score_all_genes), desc(abs_own_kME), desc(priority_gs), .by_group = TRUE) %>%
  mutate(hub_rank_all_genes = row_number()) %>%
  ungroup()

biotype_col <- pick_first_existing(hub_ranked_all, c("biotype", "biotype.x", "biotype.y", "gene_biotype", "gene_type"))

if (!is.na(biotype_col)) {
  hub_ranked_lnc <- hub_ranked_all %>%
    filter(is_lnc_biotype(.data[[biotype_col]])) %>%
    group_by(module_color) %>%
    arrange(desc(hub_score_lncRNA), desc(abs_own_kME), desc(priority_gs), .by_group = TRUE) %>%
    mutate(hub_rank_lncRNA = row_number()) %>%
    ungroup()
} else {
  hub_ranked_lnc <- tibble()
}

# Export tables

cat("Exporting ranked hub-gene tables...\n")

readr::write_csv(
  hub_ranked_all,
  file.path(output_tables_dir, "wgcna_hub_genes_ranked_focus_modules.csv")
)

top20_per_module <- hub_ranked_all %>%
  group_by(module_color) %>%
  slice_head(n = 20) %>%
  ungroup()

readr::write_csv(
  top20_per_module,
  file.path(output_tables_dir, "wgcna_top20_hub_genes_per_focus_module.csv")
)

if (nrow(hub_ranked_lnc) > 0) {
  readr::write_csv(
    hub_ranked_lnc,
    file.path(output_tables_dir, "wgcna_hub_lncRNAs_ranked_focus_modules.csv")
  )

  top10_lnc_per_module <- hub_ranked_lnc %>%
    group_by(module_color) %>%
    slice_head(n = 10) %>%
    ungroup()

  readr::write_csv(
    top10_lnc_per_module,
    file.path(output_tables_dir, "wgcna_top10_hub_lncRNAs_per_focus_module.csv")
  )
}

# Plot selected module sizes

cat("Creating module-size overview plot...\n")

plot_module_sizes <- module_sizes %>%
  mutate(
    focus_flag = ifelse(module_color %in% focus_modules, "Focus module", "Other module"),
    module_color = forcats::fct_reorder(module_color, n_genes, .desc = TRUE)
  )

p1 <- ggplot(plot_module_sizes, aes(x = module_color, y = n_genes, fill = module_color)) +
  geom_col(color = neutral_dark, linewidth = 0.3) +
  geom_point(
    data = subset(plot_module_sizes, focus_flag == "Focus module"),
    aes(x = module_color, y = n_genes),
    shape = 21,
    size = 3,
    stroke = 1,
    fill = "white",
    color = neutral_dark
  ) +
  scale_fill_identity() +
  labs(
    title = "WGCNA module sizes",
    subtitle = "Selected focus modules are marked with white points",
    x = "Module colour",
    y = "Number of genes"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(output_plots_dir, "focus_modules_module_sizes.png"),
  plot = p1,
  width = 12,
  height = 7,
  dpi = 300
)

## Plot top hub genes per module

cat("Creating top hub-gene plot...\n")

plot_top_hubs <- top20_per_module %>%
  group_by(module_color) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    label = if ("gene_symbol" %in% colnames(.)) {
      ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol)
    } else {
      gene_id
    }
  )

p2 <- ggplot(plot_top_hubs, aes(x = reorder(label, hub_score_all_genes), y = hub_score_all_genes, fill = module_color)) +
  geom_col(color = neutral_dark, width = 0.75) +
  coord_flip() +
  facet_wrap(~ module_color, scales = "free_y") +
  scale_fill_identity() +
  labs(
    title = "Top hub genes in selected WGCNA modules",
    subtitle = "Genes ranked by composite hub score within each focus module",
    x = "Gene",
    y = "Composite hub score"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_plots_dir, "top_hub_genes_focus_modules.png"),
  plot = p2,
  width = 12,
  height = 9,
  dpi = 300
)

# Plot module membership vs gene significance

# For one representative focus module, plot absolute module membership against
# absolute T2D gene significance.
best_module <- focus_modules[1]

gs_plot_col <- if ("GSabs_status_T2D" %in% colnames(hub_ranked_all)) {
  "GSabs_status_T2D"
} else if ("priority_gs" %in% colnames(hub_ranked_all)) {
  "priority_gs"
} else {
  NA_character_
}

if (!is.na(gs_plot_col)) {
  scatter_df <- hub_ranked_all %>%
    filter(module_color == best_module)

  cor_info <- safe_cor_test(scatter_df$abs_own_kME, scatter_df[[gs_plot_col]])

  p3 <- ggplot(scatter_df, aes(x = abs_own_kME, y = .data[[gs_plot_col]])) +
    geom_point(color = best_module, alpha = 0.75, size = 2.2) +
    geom_smooth(method = "lm", se = FALSE, color = neutral_dark, linewidth = 0.8) +
    labs(
      title = paste0("Module membership vs gene significance: ", best_module, " module"),
      subtitle = paste0(
        "r = ", round(cor_info["estimate"], 3),
        " | p = ", signif(cor_info["p.value"], 3),
        " (", p_to_stars(cor_info["p.value"]), ")"
      ),
      x = "Absolute module membership (|kME|)",
      y = paste0("Absolute gene significance (", gs_plot_col, ")")
    ) +
    base_theme

  ggsave(
    filename = file.path(output_plots_dir, paste0("MM_vs_GS_", best_module, ".png")),
    plot = p3,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Save R workspace

cat("Saving RData workspace...\n")

save(
  datExpr,
  traits_df,
  module_assignments,
  module_sizes,
  module_trait_summary,
  focus_modules,
  focus_module_table,
  kME_mat,
  gs_merged,
  hub_ranked_all,
  hub_ranked_lnc,
  top20_per_module,
  file = file.path(output_rdata_dir, "05C_hub_gene_prioritisation.RData")
)

# summary

cat("\n========================================\n")
cat("05C COMPLETE\n")
cat("========================================\n")
cat("Focus modules: ", paste(focus_modules, collapse = ", "), "\n", sep = "")
cat("Genes scored in focus modules: ", nrow(hub_ranked_all), "\n", sep = "")
if (nrow(hub_ranked_lnc) > 0) {
  cat("lncRNAs scored in focus modules: ", nrow(hub_ranked_lnc), "\n", sep = "")
}
cat("Outputs written to: ", output_dir, "\n", sep = "")
cat("Key outputs:\n")
cat(" - results/wgcna/hubgenes/tables/focus_modules_selected_for_hub_analysis.csv\n")
cat(" - results/wgcna/hubgenes/tables/wgcna_hub_genes_ranked_focus_modules.csv\n")
cat(" - results/wgcna/hubgenes/tables/wgcna_top20_hub_genes_per_focus_module.csv\n")
if (nrow(hub_ranked_lnc) > 0) {
  cat(" - results/wgcna/hubgenes/tables/wgcna_hub_lncRNAs_ranked_focus_modules.csv\n")
  cat(" - results/wgcna/hubgenes/tables/wgcna_top10_hub_lncRNAs_per_focus_module.csv\n")
}
cat(" - results/wgcna/hubgenes/plots/focus_modules_module_sizes.png\n")
cat(" - results/wgcna/hubgenes/plots/top_hub_genes_focus_modules.png\n")
cat(" - results/wgcna/rdata/05C_hub_gene_prioritisation.RData\n")
cat("========================================\n")
