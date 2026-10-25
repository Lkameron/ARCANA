#!/usr/bin/env Rscript
# C1_de_covariate_adjusted.R
#
# Purpose:
#   ANOVA confirmed age (p=0.0077) and BMI (p=0.0146) differ significantly
#   across the ND/IGT/T2D groups in this cohort, meaning the main DE model
#   (design = ~status, script 03) risks conflating age/BMI effects with
#   disease-status effects. This script re-fits DESeq2 with those covariates
#   included (design = ~age + bmi + status) as a fully separate, parallel
#   analysis, then compares the resulting DEG lists to the original
#   status-only model to test whether the core findings survive adjustment.
#
# Outputs:
#   results/covariate_adjusted_de/tables/*_all_genes.tsv, *_significant.tsv
#   results/covariate_adjusted_de/tables/deg_overlap_comparison.csv

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tibble)
    library(readr)
})

project_dir <- getwd()
qc_data_dir <- file.path(project_dir, "data", "processed", "qc")
out_dir <- file.path(project_dir, "results", "covariate_adjusted_de")
out_tables_dir <- file.path(out_dir, "tables")
dir.create(out_tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("C1 - Covariate-adjusted differential expression\n")
cat("========================================\n\n")

counts_file <- file.path(qc_data_dir, "counts_matrix_filtered.tsv")
metadata_file <- file.path(qc_data_dir, "sample_metadata_ND_IGT_T2D.csv")
annotation_file <- file.path(qc_data_dir, "gene_annotation_filtered.tsv")

count_matrix_filtered <- data.table::fread(counts_file) %>%
    tibble::column_to_rownames(var = colnames(.)[1]) %>%
    as.matrix()
mode(count_matrix_filtered) <- "numeric"

sample_metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)
gene_annotation <- readr::read_tsv(annotation_file, show_col_types = FALSE)

sample_metadata <- sample_metadata %>%
    dplyr::filter(sample %in% colnames(count_matrix_filtered)) %>%
    dplyr::mutate(
        status = factor(status, levels = c("ND", "IGT", "T2D")),
        age = as.numeric(age),
        bmi = as.numeric(bmi)
    )
sample_metadata <- sample_metadata[match(colnames(count_matrix_filtered), sample_metadata$sample), ]
stopifnot(identical(colnames(count_matrix_filtered), sample_metadata$sample))
stopifnot(!any(is.na(sample_metadata$age)), !any(is.na(sample_metadata$bmi)))

cat("Age range:", range(sample_metadata$age), "\n")
cat("BMI range:", round(range(sample_metadata$bmi), 1), "\n\n")

# Covariate-adjusted model: age + bmi + status
dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filtered),
    colData = as.data.frame(sample_metadata),
    design = ~age + bmi + status
)
dds$status <- relevel(dds$status, ref = "ND")

cat("Running covariate-adjusted DESeq2 (design = ~age + bmi + status)...\n")
dds <- DESeq(dds)
cat("Complete.\n\n")

annotate_results <- function(res_object, gene_annotation_table) {
    as.data.frame(res_object) %>%
        tibble::rownames_to_column(var = "gene_id") %>%
        dplyr::left_join(
            gene_annotation_table %>%
                dplyr::select(original_gene_id, gene_symbol, gene_name, biotype),
            by = c("gene_id" = "original_gene_id")
        ) %>%
        dplyr::mutate(
            gene_symbol = dplyr::if_else(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol)
        ) %>%
        dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))
}

extract_significant <- function(res_df, padj_cutoff = 0.05, lfc_cutoff = 1) {
    res_df %>%
        dplyr::filter(!is.na(padj)) %>%
        dplyr::filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
}

comparisons <- list(
    T2D_vs_ND  = c("T2D", "ND"),
    IGT_vs_ND  = c("IGT", "ND"),
    T2D_vs_IGT = c("T2D", "IGT")
)

sig_new <- list()
for (name in names(comparisons)) {
    num <- comparisons[[name]][1]
    den <- comparisons[[name]][2]
    cat("Extracting contrast:", name, "\n")

    res <- results(dds, contrast = c("status", num, den), alpha = 0.05)
    res_df <- annotate_results(res, gene_annotation)
    sig_df <- extract_significant(res_df)
    sig_new[[name]] <- sig_df

    readr::write_tsv(res_df, file.path(out_tables_dir, paste0("DESeq2_", name, "_covariate_adjusted_all_genes.tsv")))
    readr::write_tsv(sig_df, file.path(out_tables_dir, paste0("DESeq2_", name, "_covariate_adjusted_significant.tsv")))

    cat("  Significant DEGs (covariate-adjusted):", nrow(sig_df), "\n")
}

cat("\n========================================\n")
cat("Comparing to original status-only DEG lists...\n")
cat("========================================\n\n")

de_dir <- file.path(project_dir, "results", "de", "tables")
overlap_rows <- list()

for (name in names(comparisons)) {
    orig_file <- file.path(de_dir, paste0("DESeq2_", name, "_significant.tsv"))
    if (!file.exists(orig_file)) {
        cat("  [SKIP]", name, "- original file not found\n")
        next
    }
    orig_sig <- readr::read_tsv(orig_file, show_col_types = FALSE)
    new_sig <- sig_new[[name]]

    orig_genes <- unique(orig_sig$gene_id)
    new_genes  <- unique(new_sig$gene_id)

    n_orig <- length(orig_genes)
    n_new  <- length(new_genes)
    n_shared <- length(intersect(orig_genes, new_genes))
    n_lost   <- length(setdiff(orig_genes, new_genes))
    n_gained <- length(setdiff(new_genes, orig_genes))
    pct_survived <- if (n_orig > 0) round(100 * n_shared / n_orig, 1) else NA

    cat(sprintf("%-12s original=%-5d covariate-adj=%-5d shared=%-5d lost=%-5d gained=%-5d survived=%s%%\n",
                name, n_orig, n_new, n_shared, n_lost, n_gained, pct_survived))

    overlap_rows[[name]] <- data.frame(
        comparison = name,
        original_DEGs = n_orig,
        covariate_adjusted_DEGs = n_new,
        shared = n_shared,
        lost_after_adjustment = n_lost,
        gained_after_adjustment = n_gained,
        pct_original_survived = pct_survived
    )
}

overlap_df <- dplyr::bind_rows(overlap_rows)
readr::write_csv(overlap_df, file.path(out_tables_dir, "deg_overlap_comparison.csv"))

cat("\nWritten:", file.path(out_tables_dir, "deg_overlap_comparison.csv"), "\n")
cat("\nC1 COMPLETE\n")
