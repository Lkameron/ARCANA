# SCRIPT: 04_lncrna_prioritisation.R

# Package setup
required_cran <- c(
  "dplyr",
  "readr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "writexl",
  "tidyr",
  "stringr"
)

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(writexl)
library(tidyr)
library(stringr)

cat("========================================\n")
cat("lncRNA prioritisation analysis\n")
cat("========================================\n")
cat("R version: ", R.version.string, "\n", sep = "")
cat("Run date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

# Paths and output folders

project_dir <- getwd()

de_data_dir <- file.path(project_dir, "data", "processed", "de")
de_results_dir <- file.path(project_dir, "results", "de")
de_tables_dir <- file.path(de_results_dir, "tables")

lnc_results_dir <- file.path(project_dir, "results", "lncrna")
lnc_tables_dir <- file.path(lnc_results_dir, "tables")
lnc_plots_dir <- file.path(lnc_results_dir, "plots")

dir.create(lnc_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lnc_plots_dir, recursive = TRUE, showWarnings = FALSE)

# Plot settings

theme_set(
  theme_bw(base_size = 12, base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
)

comparison_colors <- c(
  "IGT vs ND" = "#345678",
  "T2D vs ND" = "#8b0000",
  "T2D vs IGT" = "#5c4033"
)

# Input files

cat("Loading input files...\n")

normalized_counts_file <- file.path(de_data_dir, "normalized_counts.tsv")

igt_vs_nd_file <- file.path(de_tables_dir, "DESeq2_IGT_vs_ND_lncRNA_only.tsv")
t2d_vs_nd_file <- file.path(de_tables_dir, "DESeq2_T2D_vs_ND_lncRNA_only.tsv")
t2d_vs_igt_file <- file.path(de_tables_dir, "DESeq2_T2D_vs_IGT_lncRNA_only.tsv")

required_files <- c(
  normalized_counts_file,
  igt_vs_nd_file,
  t2d_vs_nd_file,
  t2d_vs_igt_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

normalized_counts <- read_tsv(normalized_counts_file, show_col_types = FALSE)
lnc_igt_vs_nd <- read_tsv(igt_vs_nd_file, show_col_types = FALSE)
lnc_t2d_vs_nd <- read_tsv(t2d_vs_nd_file, show_col_types = FALSE)
lnc_t2d_vs_igt <- read_tsv(t2d_vs_igt_file, show_col_types = FALSE)

cat("Files loaded successfully.\n")
cat("  IGT vs ND lncRNAs: ", nrow(lnc_igt_vs_nd), "\n", sep = "")
cat("  T2D vs ND lncRNAs: ", nrow(lnc_t2d_vs_nd), "\n", sep = "")
cat("  T2D vs IGT lncRNAs: ", nrow(lnc_t2d_vs_igt), "\n\n", sep = "")

# Mean normalized expression

cat("Calculating mean normalized expression...\n")

stopifnot("gene_id" %in% colnames(normalized_counts))

mean_expression <- normalized_counts %>%
  mutate(
    mean_normalized_expression = rowMeans(select(., -gene_id), na.rm = TRUE)
  ) %>%
  select(gene_id, mean_normalized_expression)


rescale_01 <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(1, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

add_missing_columns <- function(df) {
  needed <- c("gene_id", "gene_symbol", "gene_name", "biotype",
              "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  
  for (nm in needed) {
    if (!nm %in% colnames(df)) {
      df[[nm]] <- NA
    }
  }
  df
}

prioritise_lncrna <- function(lnc_df, comparison_label) {
  
  lnc_df <- add_missing_columns(lnc_df)
  
  out <- lnc_df %>%
    left_join(mean_expression, by = "gene_id") %>%
    mutate(
      mean_normalized_expression = if_else(is.na(mean_normalized_expression), 0, mean_normalized_expression),
      abs_log2FC = abs(log2FoldChange),
      neg_log10_padj = if_else(!is.na(padj) & padj > 0, -log10(padj), NA_real_)
    ) %>%
    mutate(
      score_padj = rescale_01(neg_log10_padj),
      score_lfc  = rescale_01(abs_log2FC),
      score_expr = rescale_01(log1p(mean_normalized_expression))
    ) %>%
    mutate(
      prioritisation_score = (score_padj + score_lfc + score_expr) / 3
    ) %>%
    arrange(desc(prioritisation_score), desc(abs_log2FC), padj) %>%
    mutate(
      prioritisation_rank = row_number(),
      comparison = comparison_label
    ) %>%
    select(
      comparison,
      prioritisation_rank,
      gene_id,
      gene_symbol,
      gene_name,
      biotype,
      baseMean,
      log2FoldChange,
      abs_log2FC,
      lfcSE,
      stat,
      pvalue,
      padj,
      neg_log10_padj,
      mean_normalized_expression,
      score_padj,
      score_lfc,
      score_expr,
      prioritisation_score
    )
  
  return(out)
}

# Rank lncRNAs within each comparison

cat("Ranking lncRNAs within each comparison...\n")

ranked_igt_vs_nd <- prioritise_lncrna(lnc_igt_vs_nd, "IGT vs ND")
ranked_t2d_vs_nd <- prioritise_lncrna(lnc_t2d_vs_nd, "T2D vs ND")
ranked_t2d_vs_igt <- prioritise_lncrna(lnc_t2d_vs_igt, "T2D vs IGT")

all_ranked_long <- bind_rows(
  ranked_igt_vs_nd,
  ranked_t2d_vs_nd,
  ranked_t2d_vs_igt
)

#Cross-comparison summary

cat("Building cross-comparison summary...\n")

cross_comparison_summary <- all_ranked_long %>%
  mutate(
    gene_label = case_when(
      !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
      TRUE ~ gene_id
    )
  ) %>%
  group_by(gene_id, gene_symbol, gene_name, biotype, gene_label) %>%
  summarise(
    comparisons_present = n(),
    best_rank = min(prioritisation_rank, na.rm = TRUE),
    mean_rank = mean(prioritisation_rank, na.rm = TRUE),
    max_prioritisation_score = max(prioritisation_score, na.rm = TRUE),
    mean_prioritisation_score = mean(prioritisation_score, na.rm = TRUE),
    max_abs_log2FC = max(abs_log2FC, na.rm = TRUE),
    min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    mean_expression_overall = mean(mean_normalized_expression, na.rm = TRUE),
    present_in_IGT_vs_ND = any(comparison == "IGT vs ND"),
    present_in_T2D_vs_ND = any(comparison == "T2D vs ND"),
    present_in_T2D_vs_IGT = any(comparison == "T2D vs IGT"),
    .groups = "drop"
  ) %>%
  mutate(
    min_padj = ifelse(is.infinite(min_padj), NA, min_padj),
    recurrence_score = comparisons_present / 3,
    global_priority_score = (recurrence_score + max_prioritisation_score) / 2
  ) %>%
  arrange(desc(global_priority_score), best_rank, min_padj)

# T2D vs ND points of interest for network analysis

cat("Creating T2D vs ND points-of-interest table...\n")

t2d_vs_nd_points_of_interest <- ranked_t2d_vs_nd %>%
  left_join(
    cross_comparison_summary %>%
      select(
        gene_id,
        comparisons_present,
        best_rank,
        mean_rank,
        global_priority_score,
        present_in_IGT_vs_ND,
        present_in_T2D_vs_ND,
        present_in_T2D_vs_IGT
      ),
    by = "gene_id"
  ) %>%
  mutate(
    network_priority_flag = case_when(
      prioritisation_rank <= 10 & comparisons_present >= 2 ~ "high",
      prioritisation_rank <= 20 ~ "medium",
      TRUE ~ "review"
    )
  ) %>%
  arrange(prioritisation_rank)

# Write output tables

cat("Writing output tables...\n")

write_tsv(ranked_igt_vs_nd,
          file.path(lnc_tables_dir, "lncRNA_prioritised_IGT_vs_ND.tsv"))

write_tsv(ranked_t2d_vs_nd,
          file.path(lnc_tables_dir, "lncRNA_prioritised_T2D_vs_ND.tsv"))

write_tsv(ranked_t2d_vs_igt,
          file.path(lnc_tables_dir, "lncRNA_prioritised_T2D_vs_IGT.tsv"))

write_tsv(all_ranked_long,
          file.path(lnc_tables_dir, "lncRNA_prioritised_all_comparisons_long.tsv"))

write_tsv(cross_comparison_summary,
          file.path(lnc_tables_dir, "lncRNA_cross_comparison_summary.tsv"))

write_tsv(t2d_vs_nd_points_of_interest,
          file.path(lnc_tables_dir, "lncRNA_points_of_interest_T2D_vs_ND_for_networks.tsv"))

excel_output <- list(
  summary = tibble(
    comparison = c("IGT vs ND", "T2D vs ND", "T2D vs IGT"),
    n_ranked_lncrnas = c(
      nrow(ranked_igt_vs_nd),
      nrow(ranked_t2d_vs_nd),
      nrow(ranked_t2d_vs_igt)
    )
  ),
  ranked_IGT_vs_ND = ranked_igt_vs_nd,
  ranked_T2D_vs_ND = ranked_t2d_vs_nd,
  ranked_T2D_vs_IGT = ranked_t2d_vs_igt,
  cross_comparison_summary = cross_comparison_summary,
  T2D_vs_ND_points_of_interest = t2d_vs_nd_points_of_interest
)

write_xlsx(
  excel_output,
  path = file.path(project_dir, "GSE164416_lncRNA_prioritisation_summary.xlsx")
)

# Plot functions

plot_top_bar <- function(df, comparison_name, output_file, top_n = 15) {
  
  top_df <- df %>%
    slice_head(n = top_n) %>%
    mutate(
      gene_label = if_else(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id),
      gene_label = factor(gene_label, levels = rev(gene_label))
    )
  
  p <- ggplot(top_df, aes(x = gene_label, y = prioritisation_score)) +
    geom_col(fill = comparison_colors[[comparison_name]], width = 0.75) +
    coord_flip() +
    labs(
      title = paste("Top prioritised lncRNAs:", comparison_name),
      subtitle = paste("Top", top_n, "candidates ranked by significance, effect size, and expression"),
      x = NULL,
      y = "Prioritisation score"
    ) +
    theme(axis.text.y = element_text(size = 10))
  
  ggsave(output_file, p, width = 9, height = 6.5, dpi = 300)
}

plot_scatter <- function(df, comparison_name, output_file, top_labels = 10) {
  
  plot_df <- df %>%
    mutate(
      gene_label = if_else(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id),
      label_to_show = if_else(prioritisation_rank <= top_labels, gene_label, NA_character_)
    )
  
  p <- ggplot(plot_df, aes(x = mean_normalized_expression, y = neg_log10_padj)) +
    geom_point(color = comparison_colors[[comparison_name]], alpha = 0.8, size = 2.8) +
    geom_text_repel(
      aes(label = label_to_show),
      family = "serif",
      size = 3.4,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = 25
    ) +
    labs(
      title = paste("lncRNA prioritisation overview:", comparison_name),
      subtitle = "Higher significance and stronger expression support prioritisation",
      x = "Mean normalized expression",
      y = expression(-log[10]("adjusted p-value"))
    )
  
  ggsave(output_file, p, width = 8.5, height = 6.5, dpi = 300)
}

plot_recurrence <- function(df, output_file, top_n = 20) {
  
  plot_df <- df %>%
    slice_head(n = top_n) %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_label))
    )
  
  p <- ggplot(plot_df, aes(x = gene_label, y = comparisons_present, fill = global_priority_score)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_gradient(low = "#d9d9d9", high = "#8b0000") +
    labs(
      title = "Top recurrent lncRNA candidates across comparisons",
      subtitle = paste("Top", top_n, "lncRNAs ranked by recurrence and maximum prioritisation score"),
      x = NULL,
      y = "Number of comparisons",
      fill = "Global\npriority"
    )
  
  ggsave(output_file, p, width = 9, height = 7, dpi = 300)
}

# Generate plots

cat("Generating plots...\n")

plot_top_bar(
  ranked_igt_vs_nd,
  "IGT vs ND",
  file.path(lnc_plots_dir, "top_prioritised_lncRNAs_IGT_vs_ND.png"),
  top_n = 15
)

plot_top_bar(
  ranked_t2d_vs_nd,
  "T2D vs ND",
  file.path(lnc_plots_dir, "top_prioritised_lncRNAs_T2D_vs_ND.png"),
  top_n = 15
)

plot_top_bar(
  ranked_t2d_vs_igt,
  "T2D vs IGT",
  file.path(lnc_plots_dir, "top_prioritised_lncRNAs_T2D_vs_IGT.png"),
  top_n = 15
)

plot_scatter(
  ranked_t2d_vs_nd,
  "T2D vs ND",
  file.path(lnc_plots_dir, "lncRNA_prioritisation_scatter_T2D_vs_ND.png"),
  top_labels = 10
)

plot_recurrence(
  cross_comparison_summary,
  file.path(lnc_plots_dir, "lncRNA_recurrence_across_comparisons.png"),
  top_n = 20
)

# console summary

cat("\n========================================\n")
cat("lncRNA PRIORITISATION SUMMARY\n")
cat("========================================\n\n")

cat("Ranked lncRNAs by comparison:\n")
cat("  IGT vs ND: ", nrow(ranked_igt_vs_nd), "\n", sep = "")
cat("  T2D vs ND: ", nrow(ranked_t2d_vs_nd), "\n", sep = "")
cat("  T2D vs IGT: ", nrow(ranked_t2d_vs_igt), "\n\n", sep = "")

cat("Top 10 T2D vs ND points of interest for network analysis:\n")
print(
  t2d_vs_nd_points_of_interest %>%
    select(
      prioritisation_rank,
      gene_symbol,
      gene_id,
      log2FoldChange,
      padj,
      mean_normalized_expression,
      comparisons_present,
      global_priority_score,
      network_priority_flag
    ) %>%
    head(10)
)
cat("\n")

cat("Files created:\n")
cat("  - results/lncrna/tables/lncRNA_prioritised_IGT_vs_ND.tsv\n")
cat("  - results/lncrna/tables/lncRNA_prioritised_T2D_vs_ND.tsv\n")
cat("  - results/lncrna/tables/lncRNA_prioritised_T2D_vs_IGT.tsv\n")
cat("  - results/lncrna/tables/lncRNA_prioritised_all_comparisons_long.tsv\n")
cat("  - results/lncrna/tables/lncRNA_cross_comparison_summary.tsv\n")
cat("  - results/lncrna/tables/lncRNA_points_of_interest_T2D_vs_ND_for_networks.tsv\n")
cat("  - results/lncrna/plots/ (5 PNG files)\n")
cat("  - GSE164416_lncRNA_prioritisation_summary.xlsx\n")
cat("========================================\n\n")

cat("lncRNA prioritisation complete.\n")
cat("Next step: use the T2D vs ND points-of-interest table to inspect candidate lncRNAs in co-expression modules.\n")