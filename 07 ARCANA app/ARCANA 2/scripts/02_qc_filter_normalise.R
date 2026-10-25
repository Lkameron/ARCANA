# SCRIPT: Quality Control and Exploratory Analysis for GSE164416

# What this script does:
# 1. Loads count matrix, sample metadata, and gene annotation from data acquisition
# 2. Restricts the dataset to ND, IGT, and T2D samples only
# 3. Checks sample-count matrix consistency
# 4. Summarizes library sizes and sample group distribution
# 5. Filters lowly expressed genes
# 6. Creates a DESeq2 dataset for exploratory QC
# 7. Produces PCA and sample-distance heatmap plots
# 8. Exports QC tables, filtered data, and figures for downstream analysis

#  Setup 

required_packages <- c(
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tibble",
  "ggplot2",
  "pheatmap",
  "RColorBrewer"
)

bioc_packages <- c(
  "DESeq2"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)

cat("R version:", R.version.string, "\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

# Directory Structure

project_dir <- getwd()

metadata_dir <- file.path(project_dir, "data", "metadata")
processed_dir <- file.path(project_dir, "data", "processed")
qc_data_dir <- file.path(processed_dir, "qc")

results_dir <- file.path(project_dir, "results")
qc_results_dir <- file.path(results_dir, "qc")
qc_tables_dir <- file.path(qc_results_dir, "tables")
qc_plots_dir <- file.path(qc_results_dir, "plots")
logs_dir <- file.path(results_dir, "logs")

dir.create(qc_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(processed_dir, "counts_matrix.tsv")
sample_metadata_file <- file.path(metadata_dir, "sample_metadata.csv")
gene_annotation_file <- file.path(metadata_dir, "gene_annotation.tsv")

cat("QC directories created successfully.\n\n")

# Load Input Files

if (!file.exists(counts_file)) {
  stop("counts_matrix.tsv not found. Please run data_aquisition.R first.")
}

if (!file.exists(sample_metadata_file)) {
  stop("sample_metadata.csv not found. Please run data_aquisition.R first.")
}

if (!file.exists(gene_annotation_file)) {
  stop("gene_annotation.tsv not found. Please run data_aquisition.R first.")
}

cat("Loading input files...\n")

counts_df <- data.table::fread(counts_file)
sample_metadata <- readr::read_csv(sample_metadata_file, show_col_types = FALSE)
gene_annotation <- readr::read_tsv(gene_annotation_file, show_col_types = FALSE)

gene_id_col <- colnames(counts_df)[1]

count_matrix <- counts_df %>%
  tibble::column_to_rownames(var = gene_id_col) %>%
  as.matrix()

mode(count_matrix) <- "numeric"

cat("Input files loaded successfully.\n")
cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")
cat("Sample metadata rows:", nrow(sample_metadata), "\n")
cat("Gene annotation rows:", nrow(gene_annotation), "\n\n")

# Restrict to ND, IGT, and T2D Samples

cat("Subsetting to ND, IGT, and T2D samples only...\n")

sample_metadata_subset <- sample_metadata %>%
  dplyr::filter(status %in% c("ND", "IGT", "T2D")) %>%
  dplyr::mutate(
    status = factor(status, levels = c("ND", "IGT", "T2D"))
  )

retained_samples <- sample_metadata_subset$sample

count_matrix_subset <- count_matrix[, retained_samples, drop = FALSE]

if (!identical(colnames(count_matrix_subset), sample_metadata_subset$sample)) {
  stop("Sample order mismatch between count matrix and metadata after subsetting.")
}

cat("Subsetting complete.\n")
cat("Retained samples:", ncol(count_matrix_subset), "\n")
cat("ND:", sum(sample_metadata_subset$status == "ND"), "\n")
cat("IGT:", sum(sample_metadata_subset$status == "IGT"), "\n")
cat("T2D:", sum(sample_metadata_subset$status == "T2D"), "\n\n")


#  Library Size and Basic Sample QC

cat("Calculating basic QC metrics...\n")

sample_qc_table <- sample_metadata_subset %>%
  dplyr::mutate(
    library_size = colSums(count_matrix_subset)[sample],
    detected_genes = colSums(count_matrix_subset > 0)[sample],
    detected_genes_1cpm_proxy = colSums(count_matrix_subset >= 10)[sample]
  ) %>%
  dplyr::arrange(status, donor_id)

group_summary_table <- sample_qc_table %>%
  dplyr::group_by(status) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_library_size = mean(library_size),
    median_library_size = median(library_size),
    min_library_size = min(library_size),
    max_library_size = max(library_size),
    mean_detected_genes = mean(detected_genes),
    .groups = "drop"
  )

cat("Basic QC metrics calculated.\n\n")

#  Filter Lowly Expressed Genes

cat("Filtering lowly expressed genes...\n")

min_count <- 10
min_samples <- ceiling(0.20 * ncol(count_matrix_subset))

keep_genes <- rowSums(count_matrix_subset >= min_count) >= min_samples

count_matrix_filtered <- count_matrix_subset[keep_genes, , drop = FALSE]

gene_filtering_summary <- tibble::tibble(
  metric = c(
    "genes_before_filtering",
    "genes_after_filtering",
    "genes_removed",
    "min_count_threshold",
    "min_samples_threshold"
  ),
  value = c(
    nrow(count_matrix_subset),
    nrow(count_matrix_filtered),
    nrow(count_matrix_subset) - nrow(count_matrix_filtered),
    min_count,
    min_samples
  )
)

gene_annotation_filtered <- gene_annotation %>%
  dplyr::filter(original_gene_id %in% rownames(count_matrix_filtered))

cat("Filtering complete.\n")
cat("Genes before filtering:", nrow(count_matrix_subset), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")
cat("Genes removed:", nrow(count_matrix_subset) - nrow(count_matrix_filtered), "\n\n")

# SBuild DESeq2 Object for Exploratory QC

cat("Constructing DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix_filtered),
  colData = as.data.frame(sample_metadata_subset),
  design = ~ status
)

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = TRUE)

sample_metadata_subset$size_factor <- sizeFactors(dds)[sample_metadata_subset$sample]

deseq_dataset_summary <- tibble::tibble(
  metric = c(
    "samples_retained",
    "genes_retained_after_filtering",
    "design_formula"
  ),
  value = c(
    ncol(dds),
    nrow(dds),
    "~ status"
  )
)

cat("DESeq2 object created successfully.\n\n")


reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
#  PCA with Variable Contribution Barplots

cat("Generating PCA with variable contribution analysis...\n")

status_colors <- c(
  "ND" = "#012345",   # Midnight Navy
  "IGT" = "#345678",  # Deep Sky
  "T2D" = "#8b0000"   # Dark Crimson
)


biotype_colors <- c(
  "Protein-coding" = "#012345",  # Midnight Navy
  "lncRNA" = "#8b0000",          # Dark Crimson (highlighting your focus)
  "Pseudogene" = "#234567",      # Ocean Blue
  "Small RNA" = "#345678",       # Deep Sky
  "Other" = "#E6F0FA"            # Ice White
)

vsd_matrix <- assay(vsd)

pca_result <- prcomp(t(vsd_matrix), scale. = TRUE)


pca_coords <- as.data.frame(pca_result$x[, 1:9])
pca_coords$sample <- rownames(pca_coords)
pca_coords <- pca_coords %>%
  dplyr::left_join(
    sample_metadata_subset %>% dplyr::select(sample, status, donor_id),
    by = "sample"
  )

# Calculate percent variance 
percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# -------------------------------------------------------------------------
# PCA Scatter Plot (PC1 vs PC2)
# -------------------------------------------------------------------------

pca_scatter <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = status)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = status_colors) +
  labs(
    title = "PCA of Variance-Stabilized Counts",
    subtitle = "GSE164416 pancreatic islet samples (ND, IGT, T2D)",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Disease Status"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(qc_plots_dir, "pca_scatter_ND_IGT_T2D.png"),
  plot = pca_scatter,
  width = 9,
  height = 6.5,
  dpi = 300
)

# -------------------------------------------------------------------------
# Gene Loadings Analysis
# -------------------------------------------------------------------------

gene_loadings <- as.data.frame(pca_result$rotation[, 1:9])
gene_loadings$gene_id <- rownames(gene_loadings)

# Annotate with biotype and gene symbol
gene_loadings <- gene_loadings %>%
  dplyr::left_join(
    gene_annotation_filtered %>% dplyr::select(original_gene_id, biotype, gene_symbol),
    by = c("gene_id" = "original_gene_id")
  ) %>%
  dplyr::mutate(
    biotype_group = dplyr::case_when(
      biotype == "protein_coding" ~ "Protein-coding",
      biotype == "lncRNA" ~ "lncRNA",
      stringr::str_detect(biotype, "pseudogene") ~ "Pseudogene",
      biotype %in% c("miRNA", "snoRNA", "snRNA", "misc_RNA") ~ "Small RNA",
      TRUE ~ "Other"
    ),
    biotype_group = factor(
      biotype_group,
      levels = c("lncRNA", "Protein-coding", "Pseudogene", "Small RNA", "Other")
    )
  )

# long format
loadings_melted <- gene_loadings %>%
  tidyr::pivot_longer(
    cols = PC1:PC9,
    names_to = "PC",
    values_to = "loading"
  ) %>%
  dplyr::mutate(
    PC = factor(PC, levels = paste0("PC", 1:9))
  )

# -------------------------------------------------------------------------
# Top 50 Gene Loadings p
# -------------------------------------------------------------------------

top_loadings <- loadings_melted %>%
  dplyr::group_by(PC) %>%
  dplyr::slice_max(order_by = abs(loading), n = 50, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(PC, desc(abs(loading)))

top_loading_barplot <- ggplot(
  top_loadings, 
  aes(x = reorder_within(gene_symbol, loading, PC), y = loading, fill = biotype_group)
) +
  geom_col() +
  scale_x_reordered() +
  facet_wrap(~ PC, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = biotype_colors) +
  coord_flip() +
  labs(
    title = "Top 50 Gene Contributors per Principal Component",
    subtitle = "Genes with highest absolute loadings (lncRNAs highlighted in red)",
    x = "Gene Symbol",
    y = "Loading (Contribution to PC)",
    fill = "Gene Biotype"
  ) +
  theme_bw(base_size = 10, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "#E6F0FA"),
    strip.text = element_text(face = "bold", color = "#012345"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(qc_plots_dir, "pca_top50_gene_loadings.png"),
  plot = top_loading_barplot,
  width = 15,
  height = 13,
  dpi = 300
)

# -------------------------------------------------------------------------
# Export top lncRNA contributors 
# -------------------------------------------------------------------------

top_lncrna_loadings <- loadings_melted %>%
  dplyr::filter(biotype_group == "lncRNA") %>%
  dplyr::group_by(PC) %>%
  dplyr::slice_max(order_by = abs(loading), n = 20, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(PC, gene_symbol, gene_id, loading, biotype)

readr::write_csv(
  top_lncrna_loadings,
  file.path(qc_tables_dir, "top_lncRNA_contributors_per_PC.csv")
)

cat("PCA plots generated with consistent color scheme:\n")
cat("  - Scatter plot (PC1 vs PC2)\n")
cat("  - Top 50 gene loading barplot (lncRNAs in red)\n")
cat("  - Top lncRNA contributors table exported\n\n")


cat("Generating sample distance heatmap with visible labels...\n")

# ComplexHeatmap
if (!require("ComplexHeatmap")) {
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}
library(ComplexHeatmap)

if (!require("circlize")) {
  install.packages("circlize", dependencies = TRUE)
}
library(circlize)

# -------------------------------------------------------------------------
# distance matrix
# -------------------------------------------------------------------------

# Calculate sample-to-sample Euclidean distances
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

# SHORT labels: 
sample_labels <- paste0(
  sample_metadata_subset$donor_id, 
  "_", 
  sample_metadata_subset$status
)

rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

cat(sprintf("Distance matrix prepared: %d samples x %d samples\n", 
            nrow(sample_dist_matrix), ncol(sample_dist_matrix)))


custom_red_navy <- colorRampPalette(c(
  "#800020",  # Deep Burgundy (low distance = similar samples)
  "#8b0000",  # Dark Crimson
  "#850101",  # Garnet
  "#2B1B17",  # Dark Roast (transition)
  "#345678",  # Deep Sky
  "#234567",  # Ocean Blue
  "#123456",  # Abyss Blue
  "#012345"   # Midnight Navy (high distance = dissimilar samples)
))(100)

# Status colors for annotation bars 
status_colors <- c(
  "ND" = "#012345",   # Midnight Navy
  "IGT" = "#345678",  # Deep Sky
  "T2D" = "#8b0000"   # Dark Crimson
)

# -------------------------------------------------------------------------
# HEATMAP 1 (98 samples)
# -------------------------------------------------------------------------

cat("Creating full heatmap (98 samples)...\n")

# Top annotation bar (disease status colors)
column_ha <- HeatmapAnnotation(
  Status = sample_metadata_subset$status,
  col = list(Status = status_colors),
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Times", fontface = "bold"),
  annotation_legend_param = list(
    Status = list(
      title = "Disease Status",
      title_gp = gpar(fontsize = 12, fontfamily = "Times", fontface = "bold"),
      labels_gp = gpar(fontsize = 11, fontfamily = "Times")
    )
  )
)

row_ha <- rowAnnotation(
  Status = sample_metadata_subset$status,
  col = list(Status = status_colors),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "Times", fontface = "bold")
)

# heatmap object
ht_full <- Heatmap(
  matrix = sample_dist_matrix,
  name = "Distance",
  
  # color scheme
  col = custom_red_navy,
  
  # Hierarchical clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = sample_dists,
  clustering_distance_columns = sample_dists,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  # ROW 
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 5.5, fontfamily = "Times"),
  row_names_max_width = unit(3.5, "in"),  # Reserve space for labels
  
  # COLUMN 
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 5.5, fontfamily = "Times"),
  column_names_rot = 45,  # Angle labels so they don't overlap
  column_names_max_height = unit(3.5, "in"),
  
  top_annotation = column_ha,
  left_annotation = row_ha,
  
  border = TRUE,
  show_heatmap_legend = TRUE,
  
  # Legend formatting 
  heatmap_legend_param = list(
    title = "Euclidean\nDistance",
    title_gp = gpar(fontsize = 12, fontfamily = "Times", fontface = "bold"),
    labels_gp = gpar(fontsize = 10, fontfamily = "Times"),
    legend_height = unit(2, "in"),
    grid_width = unit(0.2, "in")
  ),
  
  # Title 
  column_title = "Sample-to-Sample Distance Heatmap: All 98 Samples\n(ND, IGT, T2D Pancreatic Islet RNA-seq)",
  column_title_gp = gpar(fontsize = 15, fontfamily = "Times", fontface = "bold")
)


png(
  filename = file.path(qc_plots_dir, "sample_distance_heatmap_full_98samples.png"),
  width = 5000,
  height = 4800,
  res = 300
)

draw(ht_full, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(0.5, 0.5, 0.5, 0.5), "in"))

dev.off()

cat("✓ Full heatmap saved: sample_distance_heatmap_full_98samples.png\n")
cat("  (98 samples with visible labels)\n\n")

# -------------------------------------------------------------------------
# HEATMAP 2 (36 samples)
# -------------------------------------------------------------------------

cat("Creating subset heatmap (36 samples: 12 per group)...\n")

# 12 random samples from EACH group 
set.seed(42)
n_per_group <- 12

subset_samples <- sample_metadata_subset %>%
  dplyr::group_by(status) %>%
  dplyr::slice_sample(n = n_per_group) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(status, donor_id)

cat(sprintf("Subset breakdown: %d ND, %d IGT, %d T2D = %d total\n",
            sum(subset_samples$status == "ND"),
            sum(subset_samples$status == "IGT"),
            sum(subset_samples$status == "T2D"),
            nrow(subset_samples)))

# subset of distance matrix
subset_labels <- paste0(subset_samples$donor_id, "_", subset_samples$status)
subset_dist_matrix <- sample_dist_matrix[subset_labels, subset_labels]
subset_dists <- as.dist(subset_dist_matrix)

# Subset annotations 
subset_column_ha <- HeatmapAnnotation(
  Status = subset_samples$status,
  col = list(Status = status_colors),
  annotation_name_gp = gpar(fontsize = 13, fontfamily = "Times", fontface = "bold"),
  annotation_legend_param = list(
    Status = list(
      title = "Disease Status",
      title_gp = gpar(fontsize = 14, fontfamily = "Times", fontface = "bold"),
      labels_gp = gpar(fontsize = 13, fontfamily = "Times")
    )
  )
)

subset_row_ha <- rowAnnotation(
  Status = subset_samples$status,
  col = list(Status = status_colors),
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 13, fontfamily = "Times", fontface = "bold")
)

ht_subset <- Heatmap(
  matrix = subset_dist_matrix,
  name = "Distance",
  
  col = custom_red_navy,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = subset_dists,
  clustering_distance_columns = subset_dists,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  # LARGER labels 
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontfamily = "Times"),
  row_names_max_width = unit(2.5, "in"),
  
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10, fontfamily = "Times"),
  column_names_rot = 45,
  column_names_max_height = unit(2.5, "in"),
  
  top_annotation = subset_column_ha,
  left_annotation = subset_row_ha,
  
  border = TRUE,
  
  heatmap_legend_param = list(
    title = "Euclidean\nDistance",
    title_gp = gpar(fontsize = 13, fontfamily = "Times", fontface = "bold"),
    labels_gp = gpar(fontsize = 11, fontfamily = "Times"),
    legend_height = unit(2.5, "in"),
    grid_width = unit(0.25, "in")
  ),
  
  column_title = "Sample Distance Heatmap: Subset of 36 Samples\n(12 ND + 12 IGT + 12 T2D)",
  column_title_gp = gpar(fontsize = 16, fontfamily = "Times", fontface = "bold"),
  
  # Fixed size 
  width = unit(7, "in"),
  height = unit(7, "in")
)

png(
  filename = file.path(qc_plots_dir, "sample_distance_heatmap_subset_36samples.png"),
  width = 3600,
  height = 3400,
  res = 300
)

draw(ht_subset, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(0.5, 0.5, 0.5, 0.5), "in"))

dev.off()

cat("✓ Subset heatmap saved: sample_distance_heatmap_subset_36samples.png\n")
cat("  (36 samples: 12 ND, 12 IGT, 12 T2D with clear labels)\n\n")

cat("========================================\n")
cat("HEATMAP SUMMARY\n")
cat("========================================\n")
cat("Files created:\n")
cat("1. sample_distance_heatmap_full_98samples.png\n")
cat("   → All 98 samples, labels at fontsize 5.5\n")
cat("2. sample_distance_heatmap_subset_36samples.png\n")
cat("   → 36 samples (12 per group), labels at fontsize 10\n\n")
cat("Font: Times New Roman throughout\n")
cat("Color scheme: Red (similar) to Navy (dissimilar)\n")
cat("========================================\n\n")

#  QC Plots 

cat("Generating additional QC plots...\n")

library_size_plot <- ggplot(
  sample_qc_table, 
  aes(x = status, y = library_size, fill = status)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.7, size = 2.5, color = "#2B1B17") +
  scale_fill_manual(values = status_colors) +
  labs(
    title = "Library Size by Disease Status",
    x = "Disease Status",
    y = "Total Counts per Sample"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(qc_plots_dir, "library_size_by_group.png"),
  plot = library_size_plot,
  width = 7,
  height = 5.5,
  dpi = 300
)

detected_genes_plot <- ggplot(
  sample_qc_table, 
  aes(x = status, y = detected_genes, fill = status)
) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.7, size = 2.5, color = "#2B1B17") +
  scale_fill_manual(values = status_colors) +
  labs(
    title = "Detected Genes per Sample",
    x = "Disease Status",
    y = "Number of Genes with Non-Zero Counts"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(qc_plots_dir, "detected_genes_by_group.png"),
  plot = detected_genes_plot,
  width = 7,
  height = 5.5,
  dpi = 300
)

#  QC Outputs

cat("Exporting QC outputs...\n")

filtered_counts_df <- data.frame(
  gene_id = rownames(count_matrix_filtered),
  count_matrix_filtered,
  check.names = FALSE
)

normalized_counts_df <- data.frame(
  gene_id = rownames(normalized_counts),
  normalized_counts,
  check.names = FALSE
)

readr::write_tsv(
  filtered_counts_df,
  file.path(qc_data_dir, "counts_matrix_filtered.tsv")
)

readr::write_tsv(
  data.frame(gene_id = rownames(count_matrix_subset), count_matrix_subset, check.names = FALSE),
  file.path(qc_data_dir, "counts_matrix_ND_IGT_T2D.tsv")
)

readr::write_csv(
  sample_metadata_subset,
  file.path(qc_data_dir, "sample_metadata_ND_IGT_T2D.csv")
)

readr::write_tsv(
  gene_annotation_filtered,
  file.path(qc_data_dir, "gene_annotation_filtered.tsv")
)

readr::write_csv(
  gene_filtering_summary,
  file.path(qc_tables_dir, "gene_filtering_summary.csv")
)

readr::write_csv(
  sample_qc_table,
  file.path(qc_tables_dir, "sample_qc_metrics.csv")
)

readr::write_csv(
  group_summary_table,
  file.path(qc_tables_dir, "group_qc_summary.csv")
)

readr::write_csv(
  deseq_dataset_summary,
  file.path(qc_tables_dir, "deseq_dataset_summary.csv")
)

readr::write_tsv(
  normalized_counts_df,
  file.path(qc_data_dir, "normalized_counts_vst_input.tsv")
)

saveRDS(dds, file.path(qc_data_dir, "dds_qc_ready.rds"))
saveRDS(vsd, file.path(qc_data_dir, "vsd_qc_ready.rds"))

# summary 
cat("\n========================================\n")
cat("Quality Control Complete!\n")
cat("========================================\n")
cat("Project focus: ND vs IGT vs T2D only\n")
cat("Samples retained:", ncol(count_matrix_subset), "\n")
cat("ND:", sum(sample_metadata_subset$status == "ND"), "\n")
cat("IGT:", sum(sample_metadata_subset$status == "IGT"), "\n")
cat("T2D:", sum(sample_metadata_subset$status == "T2D"), "\n")
cat("\nGenes retained after filtering:", nrow(count_matrix_filtered), "\n")
cat("Genes removed:", nrow(count_matrix_subset) - nrow(count_matrix_filtered), "\n")
cat("\nFiles written:\n")
cat("  data/processed/qc/counts_matrix_ND_IGT_T2D.tsv\n")
cat("  data/processed/qc/counts_matrix_filtered.tsv\n")
cat("  data/processed/qc/sample_metadata_ND_IGT_T2D.csv\n")
cat("  data/processed/qc/gene_annotation_filtered.tsv\n")
cat("  data/processed/qc/dds_qc_ready.rds\n")
cat("  data/processed/qc/vsd_qc_ready.rds\n")
cat("  results/qc/tables/gene_filtering_summary.csv\n")
cat("  results/qc/tables/sample_qc_metrics.csv\n")
cat("  results/qc/tables/group_qc_summary.csv\n")
cat("  results/qc/tables/deseq_dataset_summary.csv\n")
cat("  results/qc/plots/pca_plot_ND_IGT_T2D.png\n")
cat("  results/qc/plots/sample_distance_heatmap_ND_IGT_T2D.png\n")
cat("  results/qc/plots/library_size_by_group.png\n")
cat("  results/qc/plots/detected_genes_by_group.png\n")
cat("\nNext steps:\n")
cat("  1. Review PCA and heatmap for outliers\n")
cat("  2. Confirm filtering threshold is appropriate\n")
cat("  3. Proceed to differential expression analysis\n")
cat("========================================\n")

