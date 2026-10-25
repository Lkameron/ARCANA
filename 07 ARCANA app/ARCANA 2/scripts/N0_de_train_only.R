# SCRIPT: 03_differential_expression.R
# What this script does:
#   1. Run DESeq2 differential expression analysis for:
#      - IGT vs ND
#      - T2D vs ND
#      - T2D vs IGT
#   2. Export full, significant, and lncRNA-only result tables
#   3. Generate volcano plots and MA plots
#   4. Save a DESeq2 object for downstream analyses

# Package setup

required_cran <- c(
    "data.table",
    "dplyr",
    "readr",
    "tibble",
    "ggplot2",
    "ggrepel",
    "writexl"
)

required_bioc <- c(
    "DESeq2",
    "EnhancedVolcano"
)

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
library(ggplot2)
library(ggrepel)
library(writexl)
library(DESeq2)
library(EnhancedVolcano)

cat("========================================\n")
cat("Differential expression analysis\n")
cat("========================================\n")
cat("R version: ", R.version.string, "\n", sep = "")
cat("DESeq2 version: ", as.character(packageVersion("DESeq2")), "\n", sep = "")
cat("Run date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

#Paths and folders

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else "run1"
cat("RUN_ID =", RUN_ID, "\n\n")

project_dir <- getwd()

qc_data_dir <- file.path(project_dir, "data", "processed", "qc")
de_data_dir <- file.path(project_dir, "data", "processed", "de")
de_results_dir <- file.path(project_dir, "results", "nested_validation", RUN_ID, "de_train")
de_tables_dir <- file.path(de_results_dir, "tables")
de_plots_dir <- file.path(de_results_dir, "plots")

# Nested-validation leakage fix: restrict to TRAIN samples only (see
# N1_nested_split.R) so this DE fit — which feeds into hub-gene scoring via
# N3_hub_genes_train_only.R — never sees the held-out test samples' labels.
nested_split_file <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables", "train_test_split.csv")
train_sample_ids <- readr::read_csv(nested_split_file, show_col_types = FALSE) %>%
    dplyr::filter(split == "train") %>%
    dplyr::pull(sampleid)

dir.create(de_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_plots_dir, recursive = TRUE, showWarnings = FALSE)

# Plot settings

status_colors <- c(
    "ND"  = "#012345",
    "IGT" = "#345678",
    "T2D" = "#8b0000"
)

volcano_colors <- c(
    "#E6F0FA", # NS
    "#345678", # log2FC only
    "#234567", # padj only
    "#8b0000" # padj and log2FC
)

theme_set(
    theme_bw(base_size = 12, base_family = "serif") +
        theme(
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 11),
            legend.title = element_text(face = "bold"),
            panel.grid.minor = element_blank()
        )
)

# Load QC-ready inputs

cat("Loading QC-ready input files...\n")

counts_file <- file.path(qc_data_dir, "counts_matrix_filtered.tsv")
metadata_file <- file.path(qc_data_dir, "sample_metadata_ND_IGT_T2D.csv")
annotation_file <- file.path(qc_data_dir, "gene_annotation_filtered.tsv")

if (!file.exists(counts_file)) stop("Missing file: ", counts_file)
if (!file.exists(metadata_file)) stop("Missing file: ", metadata_file)
if (!file.exists(annotation_file)) stop("Missing file: ", annotation_file)

count_matrix_filtered <- data.table::fread(counts_file) %>%
    tibble::column_to_rownames(var = colnames(.)[1]) %>%
    as.matrix()

mode(count_matrix_filtered) <- "numeric"

# Restrict to TRAIN samples only (nested-validation leakage fix)
count_matrix_filtered <- count_matrix_filtered[, colnames(count_matrix_filtered) %in% train_sample_ids, drop = FALSE]
cat("Restricted to TRAIN samples only: ", ncol(count_matrix_filtered), " samples\n", sep = "")

sample_metadata <- readr::read_csv(
    metadata_file,
    show_col_types = FALSE
) %>%
    dplyr::filter(sample %in% train_sample_ids)

gene_annotation <- readr::read_tsv(
    annotation_file,
    show_col_types = FALSE
)

cat("Count matrix dimensions: ", nrow(count_matrix_filtered), " genes x ",
    ncol(count_matrix_filtered), " samples\n",
    sep = ""
)
cat("Sample group counts:\n")
print(table(sample_metadata$status))
cat("\n")

# need to prepare meta data and a DESeq2 dataset object (git hub repository - look it up)
# literature: Love et al., 2014, Genome Biology, 15:R29 [15/03/2026]

cat("Preparing metadata and DESeq2 dataset...\n")

if (!"sample" %in% colnames(sample_metadata)) {
    stop("sample_metadata must contain a column named 'sample'")
}

if (!"status" %in% colnames(sample_metadata)) {
    stop("sample_metadata must contain a column named 'status'")
}

sample_metadata <- sample_metadata %>%
    dplyr::filter(sample %in% colnames(count_matrix_filtered)) %>%
    dplyr::mutate(
        status = factor(status, levels = c("ND", "IGT", "T2D"))
    )

sample_metadata <- sample_metadata[match(colnames(count_matrix_filtered), sample_metadata$sample), ]

if (!identical(colnames(count_matrix_filtered), sample_metadata$sample)) {
    stop("Sample order mismatch between count matrix and metadata.")
}

dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filtered),
    colData = as.data.frame(sample_metadata),
    design = ~status
)

dds$status <- relevel(dds$status, ref = "ND")

cat("Reference group: ND\n")
cat("Factor levels: ", paste(levels(dds$status), collapse = ", "), "\n\n", sep = "")

# Run DESeq2

cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("DESeq2 complete.\n\n")

saveRDS(dds, file.path(de_data_dir, "dds_de_complete.rds"))


annotate_results <- function(res_object, gene_annotation_table) {
    res_df <- as.data.frame(res_object) %>%
        tibble::rownames_to_column(var = "gene_id") %>%
        dplyr::left_join(
            gene_annotation_table %>%
                dplyr::select(original_gene_id, gene_symbol, gene_name, biotype),
            by = c("gene_id" = "original_gene_id")
        ) %>%
        dplyr::mutate(
            gene_symbol = dplyr::if_else(
                is.na(gene_symbol) | gene_symbol == "",
                gene_id,
                gene_symbol
            )
        ) %>%
        dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))

    return(res_df)
}

extract_significant <- function(res_df, padj_cutoff = 0.05, lfc_cutoff = 1) {
    res_df %>%
        dplyr::filter(!is.na(padj)) %>%
        dplyr::filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
}

extract_lncrna <- function(sig_df) {
    sig_df %>%
        dplyr::filter(biotype == "lncRNA")
}

create_volcano_plot <- function(res_df, comparison_name, output_file) {
    p <- EnhancedVolcano(
        res_df,
        lab = res_df$gene_symbol,
        x = "log2FoldChange",
        y = "padj",
        title = paste("Volcano plot:", comparison_name),
        subtitle = "DESeq2 differential expression analysis",
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 2.2,
        labSize = 4.0,
        col = volcano_colors,
        colAlpha = 0.75,
        legendLabels = c("NS", "Log2FC", "padj", "padj & Log2FC"),
        legendPosition = "right",
        legendLabSize = 10,
        legendIconSize = 4,
        drawConnectors = TRUE,
        widthConnectors = 0.4,
        max.overlaps = 20,
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        border = "full",
        borderWidth = 0.8,
        borderColour = "black"
    ) +
        theme(
            text = element_text(family = "serif"),
            plot.title = element_text(size = 15, face = "bold"),
            plot.subtitle = element_text(size = 11)
        )

    ggsave(
        filename = output_file,
        plot = p,
        width = 11,
        height = 8.5,
        dpi = 300
    )
}

create_ma_plot <- function(res_object, comparison_name, output_file) {
    png(output_file, width = 2400, height = 1800, res = 300, family = "serif")
    plotMA(
        res_object,
        main = paste("MA plot:", comparison_name),
        ylim = c(-5, 5),
        cex = 0.8,
        alpha = 0.05
    )
    abline(h = c(-1, 1), col = "#8b0000", lty = 2, lwd = 1.2)
    dev.off()
}

run_comparison <- function(dds_object, numerator, denominator, gene_annotation_table) {
    res <- results(
        dds_object,
        contrast = c("status", numerator, denominator),
        alpha = 0.05
    )

    res_df <- annotate_results(res, gene_annotation_table)
    sig_df <- extract_significant(res_df)
    lncrna_df <- extract_lncrna(sig_df)

    list(
        res = res,
        full = res_df,
        sig = sig_df,
        lncrna = lncrna_df
    )
}

# Pairwise comparisons

cat("Running pairwise comparisons...\n\n")

comparison_igt_vs_nd <- run_comparison(dds, "IGT", "ND", gene_annotation)
comparison_t2d_vs_nd <- run_comparison(dds, "T2D", "ND", gene_annotation)
comparison_t2d_vs_igt <- run_comparison(dds, "T2D", "IGT", gene_annotation)

cat("IGT vs ND:\n")
cat("  Significant DE genes: ", nrow(comparison_igt_vs_nd$sig), "\n", sep = "")
cat("  Significant DE lncRNAs: ", nrow(comparison_igt_vs_nd$lncrna), "\n\n", sep = "")

cat("T2D vs ND:\n")
cat("  Significant DE genes: ", nrow(comparison_t2d_vs_nd$sig), "\n", sep = "")
cat("  Significant DE lncRNAs: ", nrow(comparison_t2d_vs_nd$lncrna), "\n\n", sep = "")

cat("T2D vs IGT:\n")
cat("  Significant DE genes: ", nrow(comparison_t2d_vs_igt$sig), "\n", sep = "")
cat("  Significant DE lncRNAs: ", nrow(comparison_t2d_vs_igt$lncrna), "\n\n", sep = "")

# Export 

cat("Writing output tables...\n")

readr::write_tsv(
    comparison_igt_vs_nd$full,
    file.path(de_tables_dir, "DESeq2_IGT_vs_ND_all_genes.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_nd$full,
    file.path(de_tables_dir, "DESeq2_T2D_vs_ND_all_genes.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_igt$full,
    file.path(de_tables_dir, "DESeq2_T2D_vs_IGT_all_genes.tsv")
)

readr::write_tsv(
    comparison_igt_vs_nd$sig,
    file.path(de_tables_dir, "DESeq2_IGT_vs_ND_significant.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_nd$sig,
    file.path(de_tables_dir, "DESeq2_T2D_vs_ND_significant.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_igt$sig,
    file.path(de_tables_dir, "DESeq2_T2D_vs_IGT_significant.tsv")
)

readr::write_tsv(
    comparison_igt_vs_nd$lncrna,
    file.path(de_tables_dir, "DESeq2_IGT_vs_ND_lncRNA_only.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_nd$lncrna,
    file.path(de_tables_dir, "DESeq2_T2D_vs_ND_lncRNA_only.tsv")
)
readr::write_tsv(
    comparison_t2d_vs_igt$lncrna,
    file.path(de_tables_dir, "DESeq2_T2D_vs_IGT_lncRNA_only.tsv")
)

summary_table <- tibble::tibble(
    comparison = c("IGT vs ND", "T2D vs ND", "T2D vs IGT"),
    total_de_genes = c(
        nrow(comparison_igt_vs_nd$sig),
        nrow(comparison_t2d_vs_nd$sig),
        nrow(comparison_t2d_vs_igt$sig)
    ),
    de_lncrnas = c(
        nrow(comparison_igt_vs_nd$lncrna),
        nrow(comparison_t2d_vs_nd$lncrna),
        nrow(comparison_t2d_vs_igt$lncrna)
    ),
    upregulated = c(
        sum(comparison_igt_vs_nd$sig$log2FoldChange > 0, na.rm = TRUE),
        sum(comparison_t2d_vs_nd$sig$log2FoldChange > 0, na.rm = TRUE),
        sum(comparison_t2d_vs_igt$sig$log2FoldChange > 0, na.rm = TRUE)
    ),
    downregulated = c(
        sum(comparison_igt_vs_nd$sig$log2FoldChange < 0, na.rm = TRUE),
        sum(comparison_t2d_vs_nd$sig$log2FoldChange < 0, na.rm = TRUE),
        sum(comparison_t2d_vs_igt$sig$log2FoldChange < 0, na.rm = TRUE)
    )
)

excel_summary <- list(
    summary = summary_table,
    IGT_vs_ND_top100 = head(comparison_igt_vs_nd$sig, 100),
    T2D_vs_ND_top100 = head(comparison_t2d_vs_nd$sig, 100),
    T2D_vs_IGT_top100 = head(comparison_t2d_vs_igt$sig, 100),
    lncRNA_IGT_vs_ND = comparison_igt_vs_nd$lncrna,
    lncRNA_T2D_vs_ND = comparison_t2d_vs_nd$lncrna,
    lncRNA_T2D_vs_IGT = comparison_t2d_vs_igt$lncrna
)

writexl::write_xlsx(
    excel_summary,
    path = file.path(project_dir, "results", "nested_validation", RUN_ID, "DE_train_only_summary.xlsx")
)

# Plots

cat("Generating volcano plots...\n")

create_volcano_plot(
    comparison_igt_vs_nd$full,
    "IGT vs ND",
    file.path(de_plots_dir, "IGT_vs_ND_volcano.png")
)

create_volcano_plot(
    comparison_t2d_vs_nd$full,
    "T2D vs ND",
    file.path(de_plots_dir, "T2D_vs_ND_volcano.png")
)

create_volcano_plot(
    comparison_t2d_vs_igt$full,
    "T2D vs IGT",
    file.path(de_plots_dir, "T2D_vs_IGT_volcano.png")
)

cat("Generating MA plots...\n")

create_ma_plot(
    comparison_igt_vs_nd$res,
    "IGT vs ND",
    file.path(de_plots_dir, "IGT_vs_ND_MA_plot.png")
)

create_ma_plot(
    comparison_t2d_vs_nd$res,
    "T2D vs ND",
    file.path(de_plots_dir, "T2D_vs_ND_MA_plot.png")
)

create_ma_plot(
    comparison_t2d_vs_igt$res,
    "T2D vs IGT",
    file.path(de_plots_dir, "T2D_vs_IGT_MA_plot.png")
)

# outputs for downstream use

normalized_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id")

readr::write_tsv(
    normalized_counts,
    file.path(de_data_dir, "normalized_counts.tsv")
)

# Final summary

cat("\n========================================\n")
cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("========================================\n\n")

cat("IGT vs ND\n")
cat("  Significant DE genes (padj < 0.05, |log2FC| > 1): ",
    nrow(comparison_igt_vs_nd$sig), "\n",
    sep = ""
)
cat("  Upregulated: ",
    sum(comparison_igt_vs_nd$sig$log2FoldChange > 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Downregulated: ",
    sum(comparison_igt_vs_nd$sig$log2FoldChange < 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Significant DE lncRNAs: ",
    nrow(comparison_igt_vs_nd$lncrna), "\n\n",
    sep = ""
)

cat("T2D vs ND\n")
cat("  Significant DE genes (padj < 0.05, |log2FC| > 1): ",
    nrow(comparison_t2d_vs_nd$sig), "\n",
    sep = ""
)
cat("  Upregulated: ",
    sum(comparison_t2d_vs_nd$sig$log2FoldChange > 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Downregulated: ",
    sum(comparison_t2d_vs_nd$sig$log2FoldChange < 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Significant DE lncRNAs: ",
    nrow(comparison_t2d_vs_nd$lncrna), "\n\n",
    sep = ""
)

cat("T2D vs IGT\n")
cat("  Significant DE genes (padj < 0.05, |log2FC| > 1): ",
    nrow(comparison_t2d_vs_igt$sig), "\n",
    sep = ""
)
cat("  Upregulated: ",
    sum(comparison_t2d_vs_igt$sig$log2FoldChange > 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Downregulated: ",
    sum(comparison_t2d_vs_igt$sig$log2FoldChange < 0, na.rm = TRUE), "\n",
    sep = ""
)
cat("  Significant DE lncRNAs: ",
    nrow(comparison_t2d_vs_igt$lncrna), "\n\n",
    sep = ""
)

cat("Files created:\n")
cat("  - data/processed/de/dds_de_complete.rds\n")
cat("  - data/processed/de/normalized_counts.tsv\n")
cat("  - results/de/tables/ (9 TSV files)\n")
cat("  - results/de/plots/ (6 PNG files)\n")
cat("  - GSE164416_differential_expression_summary.xlsx\n")
cat("========================================\n\n")

cat("Differential expression analysis complete.\n")
cat("Next step: co-expression/network analysis.\n")
