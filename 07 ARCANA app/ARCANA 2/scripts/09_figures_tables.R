################################################################################
# 09_figures_tables.R
#
# Project: Pancreatic-Islet-T2D-Coexpression
#
# Purpose:
#   Compile all publication-quality composite figures and summary tables for
#   the dissertation.  This script reads the key outputs of scripts 02–08 and
#   produces:
#     - A set of high-resolution composite figures (panels assembled with
#       cowplot / patchwork)
#     - A master Excel workbook with formatted summary sheets
#     - Standalone figure-quality PNGs and PDFs
#
# Figures produced:
#   Fig 1 — QC overview (library size, PCA, sample distances)
#   Fig 2 — Differential expression overview (volcano plots, DE gene counts)
#   Fig 3 — lncRNA prioritisation (scatter + bar charts)
#   Fig 4 — WGCNA module overview (soft-threshold, module sizes, trait heatmap)
#   Fig 5 — Hub gene spotlight (kME vs GS scatter, top hub bars)
#   Fig 6 — lncRNA network integration summary
#   Fig 7 — Enrichment summary (GO BP dot plots for focus modules)
#   Fig 8 — ML classification (ROC curves, feature importance, confusion matrix)
#
# Inputs: all prior results directories
# Outputs:
#   results/figures/  — composite PNGs and PDFs
#   GSE164416_dissertation_summary.xlsx
################################################################################

############################
# SECTION 1: Packages
############################

required_cran <- c(
  "dplyr", "readr", "tibble", "stringr", "ggplot2", "tidyr",
  "forcats", "writexl", "data.table", "scales", "patchwork",
  "ggrepel", "png", "grid", "gridExtra", "cowplot"
)

required_bioc <- c()

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(writexl)
  library(data.table)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(png)
  library(grid)
  library(gridExtra)
  library(cowplot)
})

cat("========================================\n")
cat("09 — Dissertation figures and tables\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n\n", sep = "")

############################
# SECTION 2: Paths
############################

project_dir <- getwd()

res          <- file.path(project_dir, "results")
fig_dir      <- file.path(res, "figures")
fig_dir_pdf  <- file.path(fig_dir, "pdf")
fig_dir_png  <- file.path(fig_dir, "png")

dir.create(fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir_pdf, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir_png, recursive = TRUE, showWarnings = FALSE)

############################
# SECTION 3: Global theme
############################

status_colors <- c("ND" = "#012345", "IGT" = "#345678", "T2D" = "#8B0000")

base_theme <- theme_bw(base_size = 11, base_family = "serif") +
  theme(
    text             = element_text(family = "serif"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    plot.tag         = element_text(face = "bold", size = 13)
  )

theme_set(base_theme)

############################
# SECTION 4: Helper functions
############################

safe_read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_tsv(path, show_col_types = FALSE)
}
safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}
load_png_grob <- function(path) {
  if (!file.exists(path)) return(NULL)
  img <- png::readPNG(path)
  grid::rasterGrob(img, interpolate = TRUE)
}
save_figure <- function(p, name, width = 14, height = 10, device = "both") {
  if (is.null(p)) return(invisible(NULL))
  png_path <- file.path(fig_dir_png, paste0(name, ".png"))
  pdf_path <- file.path(fig_dir_pdf, paste0(name, ".pdf"))
  if (device %in% c("png", "both")) {
    ggsave(png_path, plot = p, width = width, height = height, dpi = 300, units = "in")
    cat("  PNG:", png_path, "\n")
  }
  if (device %in% c("pdf", "both")) {
    ggsave(pdf_path, plot = p, width = width, height = height, units = "in")
    cat("  PDF:", pdf_path, "\n")
  }
}

# Wrap an external PNG as a ggplot panel
png_as_ggplot <- function(path, title = NULL) {
  grob <- load_png_grob(path)
  if (is.null(grob)) return(NULL)
  p <- ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(face = "bold", size = 11, family = "serif"))
  }
  p
}

############################
# SECTION 5: Load key data tables
############################

cat("Loading key tables...\n")

de_t2d <- safe_read_tsv(file.path(res, "de/tables/DESeq2_T2D_vs_ND_all_genes.tsv"))
de_igt <- safe_read_tsv(file.path(res, "de/tables/DESeq2_IGT_vs_ND_all_genes.tsv"))
de_t2d_igt <- safe_read_tsv(file.path(res, "de/tables/DESeq2_T2D_vs_IGT_all_genes.tsv"))

de_t2d_sig  <- safe_read_tsv(file.path(res, "de/tables/DESeq2_T2D_vs_ND_significant.tsv"))
de_igt_sig  <- safe_read_tsv(file.path(res, "de/tables/DESeq2_IGT_vs_ND_significant.tsv"))
de_t2d_igt_sig <- safe_read_tsv(file.path(res, "de/tables/DESeq2_T2D_vs_IGT_significant.tsv"))

lnc_priority <- safe_read_tsv(file.path(res, "lncrna/tables/lncRNA_prioritised_T2D_vs_ND.tsv"))
lnc_cross    <- safe_read_tsv(file.path(res, "lncrna/tables/lncRNA_cross_comparison_summary.tsv"))
lnc_network  <- safe_read_tsv(file.path(res, "lncrna/network/tables/lncRNA_module_hub_summary.tsv"))

module_assignments <- safe_read_csv(file.path(res, "wgcna/modules/wgcna_gene_module_assignments.csv"))
module_sizes       <- safe_read_csv(file.path(res, "wgcna/modules/wgcna_module_sizes.csv"))
module_trait       <- safe_read_csv(file.path(res, "wgcna/module_trait/wgcna_module_trait_correlations.csv"))
module_trait_pval  <- safe_read_csv(file.path(res, "wgcna/module_trait/wgcna_module_trait_pvalues.csv"))
hub_ranked         <- safe_read_csv(file.path(res, "wgcna/hubgenes/tables/wgcna_top20_hub_genes_per_focus_module.csv"))

qc_metrics <- safe_read_csv(file.path(res, "qc/tables/sample_qc_metrics.csv"))

############################
# SECTION 6: Figure 1 — QC overview
############################

cat("\nBuilding Figure 1 — QC overview...\n")

panel_pca     <- png_as_ggplot(file.path(res, "qc/plots/pca_scatter_ND_IGT_T2D.png"),  "A  PCA")
panel_libsize <- png_as_ggplot(file.path(res, "qc/plots/library_size_by_group.png"),   "B  Library size")
panel_genes   <- png_as_ggplot(file.path(res, "qc/plots/detected_genes_by_group.png"), "C  Detected genes")
panel_dist    <- png_as_ggplot(file.path(res, "qc/plots/sample_distance_heatmap_full_98samples.png"),
                               "D  Sample distance heatmap")

panels_qc <- Filter(Negate(is.null), list(panel_pca, panel_libsize, panel_genes, panel_dist))

if (length(panels_qc) >= 2) {
  fig1 <- wrap_plots(panels_qc, ncol = 2) +
    plot_annotation(
      title    = "Figure 1 — Quality control overview",
      subtitle = "GSE164416 pancreatic islet RNA-seq (ND, IGT, T2D)",
      theme    = theme(plot.title = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig1, "fig1_qc_overview", width = 14, height = 10)
}

############################
# SECTION 7: Figure 2 — Differential expression
############################

cat("Building Figure 2 — Differential expression...\n")

# DE gene count bar chart (built from data)
if (!is.null(de_t2d_sig) && !is.null(de_igt_sig) && !is.null(de_t2d_igt_sig)) {
  make_de_bar <- function(sig_df, full_df, label) {
    if (is.null(sig_df) || is.null(full_df)) return(NULL)
    up   <- sum(sig_df$log2FoldChange >  0, na.rm = TRUE)
    down <- sum(sig_df$log2FoldChange <  0, na.rm = TRUE)
    tibble(comparison = label, direction = c("Up", "Down"), n = c(up, down))
  }

  de_counts <- bind_rows(
    make_de_bar(de_t2d_sig,     de_t2d,     "T2D vs ND"),
    make_de_bar(de_igt_sig,     de_igt,     "IGT vs ND"),
    make_de_bar(de_t2d_igt_sig, de_t2d_igt, "T2D vs IGT")
  )

  bar_colors <- c("Up" = "#8B0000", "Down" = "#012345")

  p_de_bar <- ggplot(de_counts, aes(x = comparison, y = n, fill = direction)) +
    geom_col(position = "dodge", width = 0.65, colour = "grey20", linewidth = 0.3) +
    geom_text(aes(label = n), position = position_dodge(width = 0.65),
              vjust = -0.4, size = 3.2, family = "serif") +
    scale_fill_manual(values = bar_colors, name = "Direction") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "A  Differentially expressed genes",
      subtitle = "padj < 0.05 and |log2FC| > 1",
      x = NULL, y = "Number of DE genes"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

panel_vol_t2d <- png_as_ggplot(file.path(res, "de/plots/T2D_vs_ND_volcano.png"),   "B  T2D vs ND volcano")
panel_vol_igt <- png_as_ggplot(file.path(res, "de/plots/IGT_vs_ND_volcano.png"),   "C  IGT vs ND volcano")
panel_vol_t2i <- png_as_ggplot(file.path(res, "de/plots/T2D_vs_IGT_volcano.png"),  "D  T2D vs IGT volcano")

de_panels <- Filter(Negate(is.null), list(
  if (exists("p_de_bar")) p_de_bar else NULL,
  panel_vol_t2d, panel_vol_igt, panel_vol_t2i
))

if (length(de_panels) >= 2) {
  fig2 <- wrap_plots(de_panels, ncol = 2) +
    plot_annotation(
      title    = "Figure 2 — Differential expression analysis",
      subtitle = "DESeq2 — ND vs IGT vs T2D (GSE164416)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig2, "fig2_differential_expression", width = 14, height = 11)
}

############################
# SECTION 8: Figure 3 — lncRNA prioritisation
############################

cat("Building Figure 3 — lncRNA prioritisation...\n")

panel_lnc_scatter  <- png_as_ggplot(
  file.path(res, "lncrna/plots/lncRNA_prioritisation_scatter_T2D_vs_ND.png"),
  "A  Prioritisation scatter (T2D vs ND)"
)
panel_lnc_top_t2d  <- png_as_ggplot(
  file.path(res, "lncrna/plots/top_prioritised_lncRNAs_T2D_vs_ND.png"),
  "B  Top lncRNAs (T2D vs ND)"
)
panel_lnc_recur    <- png_as_ggplot(
  file.path(res, "lncrna/plots/lncRNA_recurrence_across_comparisons.png"),
  "C  lncRNA recurrence across comparisons"
)
panel_lnc_top_igt  <- png_as_ggplot(
  file.path(res, "lncrna/plots/top_prioritised_lncRNAs_IGT_vs_ND.png"),
  "D  Top lncRNAs (IGT vs ND)"
)

lnc_panels <- Filter(Negate(is.null),
  list(panel_lnc_scatter, panel_lnc_top_t2d, panel_lnc_recur, panel_lnc_top_igt))

if (length(lnc_panels) >= 2) {
  fig3 <- wrap_plots(lnc_panels, ncol = 2) +
    plot_annotation(
      title    = "Figure 3 — lncRNA prioritisation",
      subtitle = "Differentially expressed lncRNAs across disease-status comparisons",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig3, "fig3_lncrna_prioritisation", width = 14, height = 10)
}

############################
# SECTION 9: Figure 4 — WGCNA overview
############################

cat("Building Figure 4 — WGCNA overview...\n")

panel_threshold <- png_as_ggplot(
  file.path(res, "wgcna/plots/wgcna_soft_threshold_scale_free_fit.png"),
  "A  Soft-threshold selection"
)
panel_mod_sizes <- png_as_ggplot(
  file.path(res, "wgcna/modules/wgcna_module_sizes_barplot.png"),
  "B  Module sizes"
)
panel_dendrogram <- png_as_ggplot(
  file.path(res, "wgcna/modules/wgcna_gene_dendrogram_with_module_colors.png"),
  "C  Gene dendrogram"
)
panel_trait_heat <- png_as_ggplot(
  file.path(res, "wgcna/module_trait/wgcna_module_trait_heatmap.png"),
  "D  Module–trait heatmap"
)

wgcna_panels <- Filter(Negate(is.null),
  list(panel_threshold, panel_mod_sizes, panel_dendrogram, panel_trait_heat))

if (length(wgcna_panels) >= 2) {
  fig4 <- wrap_plots(wgcna_panels, ncol = 2) +
    plot_annotation(
      title    = "Figure 4 — WGCNA co-expression network",
      subtitle = "Weighted gene co-expression network analysis — focus on disease-associated modules",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig4, "fig4_wgcna_overview", width = 14, height = 11)
}

############################
# SECTION 10: Figure 5 — Hub gene spotlight
############################

cat("Building Figure 5 — Hub gene spotlight...\n")

if (!is.null(hub_ranked)) {
  focus_modules <- unique(hub_ranked$module_color)

  # Top 10 hub genes per module — lollipop chart
  plot_top_hubs <- hub_ranked %>%
    group_by(module_color) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      label = if ("gene_symbol.x" %in% colnames(.)) {
        ifelse(is.na(gene_symbol.x) | gene_symbol.x == "" | gene_symbol.x == gene_id,
               gene_id, gene_symbol.x)
      } else gene_id,
      label = reorder(label, hub_score_all_genes)
    )

  p_hub_lollipop <- ggplot(plot_top_hubs,
                            aes(x = hub_score_all_genes, y = label, colour = module_color)) +
    geom_segment(aes(x = 0, xend = hub_score_all_genes, y = label, yend = label),
                 linewidth = 0.6) +
    geom_point(size = 2.8) +
    scale_colour_identity() +
    facet_wrap(~ module_color, scales = "free_y", ncol = 2) +
    labs(
      title    = "A  Top hub genes per focus module",
      subtitle = "Ranked by composite hub score (kME + GS + DE evidence)",
      x        = "Composite hub score",
      y        = NULL
    ) +
    base_theme +
    theme(
      legend.position = "none",
      strip.text      = element_text(face = "bold"),
      axis.text.y     = element_text(size = 7.5)
    )

  panel_mm_gs <- png_as_ggplot(
    file.path(res, "wgcna/hubgenes/plots/MM_vs_GS_brown.png"),
    "B  Module membership vs gene significance (brown)"
  )

  hub_panels <- Filter(Negate(is.null), list(p_hub_lollipop, panel_mm_gs))

  if (length(hub_panels) >= 1) {
    fig5 <- if (length(hub_panels) == 2) {
      hub_panels[[1]] + hub_panels[[2]] +
        plot_layout(ncol = 2, widths = c(2, 1))
    } else {
      hub_panels[[1]]
    }
    fig5 <- fig5 +
      plot_annotation(
        title    = "Figure 5 — Hub gene prioritisation",
        subtitle = "Top-ranked hub genes and module membership characteristics",
        theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                         plot.subtitle = element_text(size = 10, family = "serif"))
      )
    save_figure(fig5, "fig5_hub_genes", width = 14, height = 10)
  }
}

############################
# SECTION 11: Figure 6 — lncRNA network integration
############################

cat("Building Figure 6 — lncRNA network integration...\n")

if (!is.null(lnc_network)) {
  # Bar chart: lncRNAs by module with highlight on top candidates
  lnc_net_plot <- lnc_network %>%
    arrange(desc(max_prioritisationscore)) %>%
    slice_head(n = 25) %>%
    mutate(
      label = if ("gene_symbol.x" %in% colnames(.)) {
        ifelse(is.na(gene_symbol.x) | gene_symbol.x == "", gene_id, gene_symbol.x)
      } else gene_id,
      module_present = as.character(module_color)
    )

  # module_color column from lnc_network may come from the integration table
  mc_col <- intersect(c("module_color", "module"), colnames(lnc_net_plot))[1]

  p_lnc_net <- ggplot(lnc_net_plot,
    aes(x = reorder(label, max_prioritisationscore),
        y = max_prioritisationscore,
        fill = .data[[mc_col]])) +
    geom_col(colour = "grey25", linewidth = 0.25) +
    scale_fill_identity() +
    coord_flip() +
    labs(
      title    = "A  Top lncRNA candidates by prioritisation score",
      subtitle = "Coloured by WGCNA module membership",
      x        = NULL,
      y        = "Max prioritisation score"
    ) +
    base_theme +
    theme(legend.position = "none", axis.text.y = element_text(size = 8))

  fig6 <- p_lnc_net +
    plot_annotation(
      title    = "Figure 6 — lncRNA network integration",
      subtitle = "lncRNA candidates prioritised across DE comparisons and WGCNA modules",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig6, "fig6_lncrna_network_integration", width = 10, height = 8)
}

############################
# SECTION 12: Figure 7 — Enrichment summary
############################

cat("Building Figure 7 — Enrichment summary...\n")

enrich_cross <- file.path(res, "enrichment/plots/GO_BP_cross_module_summary.png")
enrich_de    <- file.path(res, "enrichment/plots/GO_BP_DE_comparisons_summary.png")

panel_enrich_cross <- png_as_ggplot(enrich_cross, "A  GO BP — focus modules")
panel_enrich_de    <- png_as_ggplot(enrich_de,    "B  GO BP — DE comparisons")

enrich_panels <- Filter(Negate(is.null), list(panel_enrich_cross, panel_enrich_de))

if (length(enrich_panels) >= 1) {
  fig7 <- wrap_plots(enrich_panels, ncol = 1) +
    plot_annotation(
      title    = "Figure 7 — Functional enrichment analysis",
      subtitle = "GO Biological Process enrichment for WGCNA focus modules and DE comparisons",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig7, "fig7_enrichment_summary", width = 13, height = 14)
}

############################
# SECTION 13: Figure 8 — Machine learning
############################

cat("Building Figure 8 — ML classification...\n")

# The hub gene expression heatmap (source PNG 4169x9270, extreme portrait
# aspect ratio) is deliberately excluded from this composite and saved as its
# own standalone figure below (Figure 9) — forcing it into a uniform grid
# alongside near-square panels squashed every panel to match its height.
panel_roc     <- png_as_ggplot(file.path(res, "ml/plots/binary_ROC_curves.png"),       "A  ROC curves")
panel_cm      <- png_as_ggplot(file.path(res, "ml/plots/confusion_matrix_RandomForest.png"),
                               "B  Confusion matrix")
panel_pca_ml  <- png_as_ggplot(file.path(res, "ml/plots/hub_gene_PCA_biplot.png"),     "C  PCA biplot")
panel_imp     <- png_as_ggplot(file.path(res, "ml/plots/RF_feature_importances.png"),  "D  Feature importances")
panel_perm    <- png_as_ggplot(file.path(res, "ml/plots/permutation_test.png"),        "E  Permutation test")

# Group panels by native aspect ratio so each row is visually balanced:
# row 1 = three near-square panels; row 2 = two wide/short panels.
top_row    <- Filter(Negate(is.null), list(panel_roc, panel_cm, panel_pca_ml))
bottom_row <- Filter(Negate(is.null), list(panel_imp, panel_perm))

row_plots <- list()
if (length(top_row) > 0)    row_plots$top    <- wrap_plots(top_row, nrow = 1)
if (length(bottom_row) > 0) row_plots$bottom <- wrap_plots(bottom_row, nrow = 1)

if (length(row_plots) >= 1) {
  fig8 <- wrap_plots(row_plots, ncol = 1, heights = rep(1, length(row_plots))) +
    plot_annotation(
      title    = "Figure 8 — Machine-learning classification of T2D",
      subtitle = "Random Forest, Logistic Regression, SVM — hub gene features (5-fold CV)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig8, "fig8_ml_classification", width = 14, height = 11)
}

############################
# SECTION 13B: Figure 9 — Hub gene expression heatmap (standalone)
############################

cat("Building Figure 9 — Hub gene expression heatmap (standalone)...\n")

panel_heatmap_standalone <- png_as_ggplot(
  file.path(res, "ml/plots/hub_gene_expression_heatmap.png"),
  NULL
)

if (!is.null(panel_heatmap_standalone)) {
  fig9 <- panel_heatmap_standalone +
    plot_annotation(
      title    = "Figure 9 — Hub gene expression across samples",
      subtitle = "Z-score-normalised expression, samples sorted by disease status (ND / IGT / T2D)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14, family = "serif"),
                       plot.subtitle = element_text(size = 10, family = "serif"))
    )
  save_figure(fig9, "fig9_hub_gene_expression_heatmap", width = 10, height = 16)
}

############################
# SECTION 14: Summary tables — master Excel workbook
############################

cat("\nCompiling master Excel workbook...\n")

sheets <- list()

# QC summary
if (!is.null(qc_metrics)) {
  sheets[["QC_sample_metrics"]] <- qc_metrics
}

# DE tables
de_summary_list <- list(
  T2D_vs_ND_sig  = de_t2d_sig,
  IGT_vs_ND_sig  = de_igt_sig,
  T2D_vs_IGT_sig = de_t2d_igt_sig
)
for (nm in names(de_summary_list)) {
  df <- de_summary_list[[nm]]
  if (!is.null(df) && nrow(df) > 0) {
    sheets[[paste0("DE_", nm)]] <- df %>%
      arrange(padj) %>%
      select(any_of(c("gene_id", "gene_symbol", "gene_name", "biotype",
                      "log2FoldChange", "padj", "pvalue", "baseMean")))
  }
}

# lncRNA prioritisation
if (!is.null(lnc_priority) && nrow(lnc_priority) > 0) {
  sheets[["lncRNA_T2D_vs_ND"]] <- lnc_priority %>% arrange(desc(prioritisation_score))
}
if (!is.null(lnc_cross) && nrow(lnc_cross) > 0) {
  sheets[["lncRNA_cross_comparison"]] <- lnc_cross
}
if (!is.null(lnc_network) && nrow(lnc_network) > 0) {
  sheets[["lncRNA_network_summary"]] <- lnc_network
}

# WGCNA module sizes
if (!is.null(module_sizes)) {
  sheets[["WGCNA_module_sizes"]] <- module_sizes
}

# Hub genes
if (!is.null(hub_ranked)) {
  sheets[["Hub_genes_top20"]] <- hub_ranked %>%
    arrange(module_color, hub_rank_all_genes) %>%
    select(any_of(c("module_color", "gene_id", "gene_symbol", "gene_name", "biotype",
                    "hub_rank_all_genes", "abs_own_kME", "priority_gs", "hub_score_all_genes")))
}

# ML metrics
ml_binary_path <- file.path(res, "ml/tables/binary_cv_metrics.csv")
ml_multi_path  <- file.path(res, "ml/tables/multiclass_cv_metrics.csv")
ml_imp_path    <- file.path(res, "ml/tables/RF_permutation_importances.csv")

if (file.exists(ml_binary_path)) sheets[["ML_binary_metrics"]]  <- readr::read_csv(ml_binary_path, show_col_types = FALSE)
if (file.exists(ml_multi_path))  sheets[["ML_multiclass"]]       <- readr::read_csv(ml_multi_path,  show_col_types = FALSE)
if (file.exists(ml_imp_path))    sheets[["ML_feature_importance"]] <- readr::read_csv(ml_imp_path,  show_col_types = FALSE)

# Write workbook
if (length(sheets) > 0) {
  writexl::write_xlsx(sheets,
    file.path(project_dir, "GSE164416_dissertation_summary.xlsx"))
  cat("  Master Excel workbook written: GSE164416_dissertation_summary.xlsx\n")
}

############################
# SECTION 15: Supplementary table — all DE results combined
############################

cat("Writing supplementary DE table...\n")

de_all_combined <- bind_rows(
  if (!is.null(de_t2d))     de_t2d     %>% mutate(comparison = "T2D_vs_ND")  else NULL,
  if (!is.null(de_igt))     de_igt     %>% mutate(comparison = "IGT_vs_ND")  else NULL,
  if (!is.null(de_t2d_igt)) de_t2d_igt %>% mutate(comparison = "T2D_vs_IGT") else NULL
)

if (nrow(de_all_combined) > 0) {
  readr::write_tsv(
    de_all_combined,
    file.path(res, "figures", "supplementary_all_DE_genes.tsv")
  )
}

############################
# SECTION 16: Final summary
############################

cat("\n========================================\n")
cat("09 COMPLETE\n")
cat("========================================\n")
cat("Figures written to: results/figures/\n")
cat("  png/ — high-resolution PNGs\n")
cat("  pdf/ — vector PDFs\n")
cat("Master Excel: GSE164416_dissertation_summary.xlsx\n")
cat("\nFigures produced:\n")
for (fig in c("fig1_qc_overview", "fig2_differential_expression",
              "fig3_lncrna_prioritisation", "fig4_wgcna_overview",
              "fig5_hub_genes", "fig6_lncrna_network_integration",
              "fig7_enrichment_summary", "fig8_ml_classification")) {
  png_p <- file.path(fig_dir_png, paste0(fig, ".png"))
  status <- if (file.exists(png_p)) "OK" else "MISSING (check upstream outputs)"
  cat(sprintf("  %-42s %s\n", fig, status))
}
cat("========================================\n")
