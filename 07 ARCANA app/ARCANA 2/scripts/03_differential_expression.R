# SCRIPT: 03_differential_expression.R (ARCANA — generalised)
#
# Runs DESeq2 differential expression for all comparisons defined in
# config/project.yaml. No hardcoded group names — reads everything from
# environment variables set by the backend (via config_utils.R).
#
# Outputs per comparison:
#   results/de/tables/DESeq2_<name>_all_genes.tsv
#   results/de/tables/DESeq2_<name>_significant.tsv
#   results/de/tables/DESeq2_<name>_lncRNA_only.tsv  (if lncrna_enabled)
#   results/de/plots/<name>_volcano.png
#   results/de/plots/<name>_MA_plot.png

source(file.path(Sys.getenv("ARCANA_SCRIPTS_DIR", "scripts"), "config_utils.R"))

required_cran <- c("data.table", "dplyr", "readr", "tibble", "ggplot2", "ggrepel", "writexl")
required_bioc <- c("DESeq2", "EnhancedVolcano")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_cran) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
for (pkg in required_bioc) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tibble)
  library(ggplot2); library(ggrepel); library(DESeq2); library(EnhancedVolcano)
})

cat("========================================\n")
cat("ARCANA — Differential Expression\n")
cat("Project:", arcana_cfg$project_name, "\n")
cat("========================================\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
cat("Run date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

# ── Paths ─────────────────────────────────────────────────────────────────────

qc_data_dir  <- file.path("data", "processed", "qc")
de_data_dir  <- file.path("data", "processed", "de")
de_tables    <- file.path(arcana_cfg$de_dir, "tables")
de_plots     <- file.path(arcana_cfg$de_dir, "plots")

dir.create(de_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_tables,   recursive = TRUE, showWarnings = FALSE)
dir.create(de_plots,    recursive = TRUE, showWarnings = FALSE)

# ── Load inputs ───────────────────────────────────────────────────────────────

cat("Loading QC-ready input files...\n")

counts_file     <- file.path(qc_data_dir, "counts_matrix_filtered.tsv")
annotation_file <- file.path(qc_data_dir, "gene_annotation_filtered.tsv")

# Metadata: prefer the QC-output version, fall back to the raw config path
metadata_candidates <- c(
  list.files(qc_data_dir, pattern = "sample_metadata.*\\.csv", full.names = TRUE),
  arcana_cfg$metadata_path
)
metadata_file <- metadata_candidates[file.exists(metadata_candidates)][1]

if (!file.exists(counts_file))   stop("Missing: ", counts_file)
if (!file.exists(metadata_file)) stop("Missing metadata file (checked: ", paste(metadata_candidates, collapse = ", "), ")")
if (!file.exists(annotation_file)) stop("Missing: ", annotation_file)

count_matrix <- data.table::fread(counts_file) %>%
  tibble::column_to_rownames(var = colnames(.)[1]) %>%
  as.matrix()
mode(count_matrix) <- "numeric"

sample_metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)
gene_annotation <- readr::read_tsv(annotation_file, show_col_types = FALSE)

cat("Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("Condition distribution:\n")
print(table(sample_metadata[[arcana_cfg$condition_col]]))
cat("\n")

# ── Build DESeq2 dataset ───────────────────────────────────────────────────────

cat("Preparing DESeq2 dataset...\n")

sid <- arcana_cfg$sample_id_col
cond <- arcana_cfg$condition_col

if (!sid %in% colnames(sample_metadata))  stop("sample_metadata missing column: ", sid)
if (!cond %in% colnames(sample_metadata)) stop("sample_metadata missing column: ", cond)

# Collect all groups mentioned in comparisons so we can set factor levels
all_groups <- unique(unlist(lapply(arcana_cfg$comparisons, function(c) c(c$group_b, c$group_a))))

sample_metadata <- sample_metadata %>%
  dplyr::filter(.data[[sid]] %in% colnames(count_matrix)) %>%
  dplyr::mutate(!!cond := factor(.data[[cond]], levels = all_groups))

sample_metadata <- sample_metadata[match(colnames(count_matrix), sample_metadata[[sid]]), ]

if (!identical(colnames(count_matrix), sample_metadata[[sid]])) {
  stop("Sample order mismatch between count matrix and metadata.")
}

# Build design formula (condition + optional covariates)
covariates <- arcana_cfg$de_covariates
covariates <- covariates[covariates %in% colnames(sample_metadata)]
design_vars <- c(covariates, cond)
design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))
cat("DESeq2 design:", deparse(design_formula), "\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = as.data.frame(sample_metadata),
  design    = design_formula
)

# Reference = first group_b in the comparison list
ref_group <- arcana_cfg$comparisons[[1]]$group_b
dds[[cond]] <- relevel(dds[[cond]], ref = ref_group)
cat("Reference group:", ref_group, "\n")
cat("Factor levels:", paste(levels(dds[[cond]]), collapse = " > "), "\n\n")

cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("DESeq2 complete.\n\n")

saveRDS(dds, file.path(de_data_dir, "dds_de_complete.rds"))

# ── Helper functions ──────────────────────────────────────────────────────────

annotate_results <- function(res_object, gene_annotation_table) {
  as.data.frame(res_object) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(
      gene_annotation_table %>% dplyr::select(original_gene_id, gene_symbol, gene_name, biotype),
      by = c("gene_id" = "original_gene_id")
    ) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol)) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
}

extract_significant <- function(res_df) {
  res_df %>%
    dplyr::filter(!is.na(padj), padj < arcana_cfg$padj, abs(log2FoldChange) > arcana_cfg$log2fc)
}

extract_lncrna <- function(sig_df) {
  sig_df %>% dplyr::filter(biotype %in% arcana_cfg$lncrna_biotypes)
}

run_comparison <- function(dds_obj, comp) {
  cat("  Running:", comp$name, "(", comp$group_a, "vs", comp$group_b, ")\n")
  res <- results(dds_obj, contrast = c(cond, comp$group_a, comp$group_b), alpha = arcana_cfg$padj)
  res_df  <- annotate_results(res, gene_annotation)
  sig_df  <- extract_significant(res_df)
  lnc_df  <- extract_lncrna(sig_df)
  cat("    Significant genes:", nrow(sig_df),
      "| lncRNA:", nrow(lnc_df), "\n")
  list(res = res, all = res_df, sig = sig_df, lncrna = lnc_df)
}

save_volcano <- function(res_df, comp_label, output_file) {
  p <- EnhancedVolcano(
    res_df,
    lab      = res_df$gene_symbol,
    x        = "log2FoldChange",
    y        = "padj",
    title    = paste("Volcano:", comp_label),
    subtitle = paste0("DESeq2 | padj<", arcana_cfg$padj, " | |log2FC|>", arcana_cfg$log2fc),
    pCutoff  = arcana_cfg$padj,
    FCcutoff = arcana_cfg$log2fc,
    pointSize = 2.2,
    labSize   = 4.0,
    drawConnectors = TRUE,
    max.overlaps   = 20,
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  )
  ggsave(output_file, plot = p, width = 11, height = 8.5, dpi = arcana_cfg$figures_dpi)
}

save_ma <- function(res_obj, comp_label, output_file) {
  png(output_file, width = 2400, height = 1800, res = arcana_cfg$figures_dpi)
  plotMA(res_obj, main = paste("MA plot:", comp_label), ylim = c(-5, 5), alpha = arcana_cfg$padj)
  abline(h = c(-1, 1), col = "#8b0000", lty = 2, lwd = 1.2)
  dev.off()
}

# ── Run all comparisons ───────────────────────────────────────────────────────

cat("Running", length(arcana_cfg$comparisons), "comparisons...\n\n")

results_list <- lapply(arcana_cfg$comparisons, function(comp) {
  r <- run_comparison(dds, comp)
  name <- comp$name

  # Save TSV tables
  readr::write_tsv(r$all,    file.path(de_tables, paste0("DESeq2_", name, "_all_genes.tsv")))
  readr::write_tsv(r$sig,    file.path(de_tables, paste0("DESeq2_", name, "_significant.tsv")))

  if (arcana_cfg$lncrna_enabled) {
    readr::write_tsv(r$lncrna, file.path(de_tables, paste0("DESeq2_", name, "_lncRNA_only.tsv")))
  }

  # Save plots
  tryCatch(
    save_volcano(r$all, comp$label %||% name, file.path(de_plots, paste0(name, "_volcano.png"))),
    error = function(e) message("  [WARNING] Volcano plot failed: ", e$message)
  )
  tryCatch(
    save_ma(r$res, comp$label %||% name, file.path(de_plots, paste0(name, "_MA_plot.png"))),
    error = function(e) message("  [WARNING] MA plot failed: ", e$message)
  )

  r
})
names(results_list) <- sapply(arcana_cfg$comparisons, `[[`, "name")

# ── Summary table ─────────────────────────────────────────────────────────────

summary_df <- data.frame(
  comparison   = arcana_cfg$comparison_names,
  n_all        = sapply(results_list, function(r) nrow(r$all)),
  n_significant = sapply(results_list, function(r) nrow(r$sig)),
  n_lncrna     = sapply(results_list, function(r) nrow(r$lncrna)),
  stringsAsFactors = FALSE
)

cat("\n── Summary ──────────────────────────────────────\n")
print(summary_df)

readr::write_csv(summary_df, file.path(de_tables, "de_summary.csv"))

# ── Excel export ──────────────────────────────────────────────────────────────

excel_sheets <- c(
  lapply(seq_along(results_list), function(i) {
    list(
      paste0(names(results_list)[i], "_top100") = head(results_list[[i]]$sig, 100)
    )
  })
)

wb_list <- setNames(
  lapply(results_list, function(r) head(r$sig, 100)),
  paste0(names(results_list), "_top100")
)
writexl::write_xlsx(wb_list, file.path(de_tables, "de_significant_top100.xlsx"))

cat("\nDifferential expression complete.\n")
cat("Tables saved to:", de_tables, "\n")
cat("Plots saved to: ", de_plots, "\n\n")
