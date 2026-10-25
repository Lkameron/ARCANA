# SCRIPT: 05B_wgcna_construct_network.R

# Package setup
required_cran <- c(
  "data.table",
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "ggplot2"
)

required_bioc <- c(
  "WGCNA",
  "ComplexHeatmap",
  "circlize"
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
library(stringr)
library(ggplot2)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()
allowWGCNAThreads()

cat("========================================\n")
cat("05B - Construct WGCNA network\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

# Define paths

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else "run1"
cat("RUN_ID =", RUN_ID, "\n\n")

project_dir <- getwd()

input_expr_file   <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables", "wgcna_expression_input_TRAIN.tsv")
input_traits_file <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables", "wgcna_traits_input_TRAIN.csv")
input_gene_annot  <- file.path(project_dir, "data", "metadata", "gene_annotation.tsv")

output_dir         <- file.path(project_dir, "results", "nested_validation", RUN_ID, "wgcna_train")
output_tables_dir  <- file.path(output_dir, "tables")
output_plots_dir   <- file.path(output_dir, "plots")
output_rdata_dir   <- file.path(output_dir, "rdata")
output_modules_dir <- file.path(output_dir, "modules")
output_trait_dir   <- file.path(output_dir, "module_trait")

for (dir_path in c(output_tables_dir, output_plots_dir, output_rdata_dir, output_modules_dir, output_trait_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# plotting style: for aesthestics 

status_colors <- c(
  "ND"  = "#012345",
  "IGT" = "#345678",
  "T2D" = "#8b0000"
)

neutral_dark  <- "#2B1B17"

base_theme_times <- theme_bw(base_size = 12, base_family = "serif") +
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

check_file_exists_nonempty <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) {
    stop("File exists but is empty: ", path)
  }
}

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

safe_cor_p <- function(x, y) {
  out <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
  if (is.null(out)) return(c(NA_real_, NA_real_))
  c(unname(out$estimate), out$p.value)
}

make_numeric_trait <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

# Load input data

for (f in c(input_expr_file, input_traits_file)) {
  check_file_exists_nonempty(f)
}

cat("Reading WGCNA inputs...\n")
expr_df <- data.table::fread(input_expr_file, data.table = FALSE)
traits_df <- readr::read_csv(input_traits_file, show_col_types = FALSE)

cat("Expression input dimensions: ", nrow(expr_df), " rows x ", ncol(expr_df), " columns\n", sep = "")
cat("Trait input dimensions: ", nrow(traits_df), " rows x ", ncol(traits_df), " columns\n\n", sep = "")

if (!("sampleid" %in% colnames(expr_df))) {
  stop("wgcna_expression_input.tsv must contain a 'sampleid' column.")
}
if (!("sampleid" %in% colnames(traits_df))) {
  stop("wgcna_traits_input.csv must contain a 'sampleid' column.")
}

# Prepare expression matrix

cat("Preparing expression matrix...\n")

datExpr <- expr_df %>%
  tibble::column_to_rownames("sampleid") %>%
  as.data.frame()

datExpr[] <- lapply(datExpr, as.numeric)
traits_df$sampleid <- as.character(traits_df$sampleid)
rownames(datExpr) <- as.character(rownames(datExpr))

common_samples <- intersect(rownames(datExpr), traits_df$sampleid)
if (length(common_samples) < 10) {
  stop("Too few overlapping samples between expression and trait inputs.")
}

datExpr <- datExpr[common_samples, , drop = FALSE]
datTraits <- traits_df %>%
  filter(sampleid %in% common_samples) %>%
  dplyr::slice(match(common_samples, sampleid))

stopifnot(all(rownames(datExpr) == datTraits$sampleid))

cat("Matched WGCNA matrix: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes\n\n", sep = "")

# WGCNA QC

cat("Running final goodSamplesGenes() check...\n")
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  datTraits <- datTraits %>%
    filter(sampleid %in% rownames(datExpr)) %>%
    dplyr::slice(match(rownames(datExpr), sampleid))
  cat("After QC filtering: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes\n\n", sep = "")
} else {
  cat("All samples and genes passed final QC.\n\n")
}

# Soft-threshold power selection

cat("Selecting soft-thresholding power...\n")

powers <- c(1:10, seq(12, 30, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1),
  verbose = 5
)

fit_df <- tibble(
  Power = sft$fitIndices[, 1],
  SFT_R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  Slope = sft$fitIndices[, 3],
  MeanConnectivity = sft$fitIndices[, 5],
  MedianConnectivity = sft$fitIndices[, 6],
  MaxConnectivity = sft$fitIndices[, 7]
)

candidate_powers <- fit_df %>%
  filter(SFT_R2 >= 0.80)

if (nrow(candidate_powers) > 0) {
  soft_power <- min(candidate_powers$Power)
} else {
  soft_power <- fit_df$Power[which.max(fit_df$SFT_R2)]
}

cat("Chosen soft-thresholding power: ", soft_power, "\n\n", sep = "")

readr::write_csv(fit_df, file.path(output_tables_dir, "wgcna_soft_threshold_results.csv"))

p_sft_fit <- ggplot(fit_df, aes(x = Power, y = SFT_R2)) +
  geom_point(size = 3, color = status_colors["ND"]) +
  geom_line(linewidth = 0.7, color = status_colors["ND"]) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = status_colors["T2D"]) +
  geom_vline(xintercept = soft_power, linetype = "dashed", color = status_colors["IGT"]) +
  geom_text(aes(label = Power), nudge_y = 0.03, size = 3.2, family = "serif") +
  labs(
    title = "Scale-free topology model fit",
    subtitle = paste0("Chosen soft-thresholding power = ", soft_power),
    x = "Soft-thresholding power",
    y = "Signed R^2"
  ) +
  coord_cartesian(ylim = c(min(0, min(fit_df$SFT_R2, na.rm = TRUE)), 1)) +
  base_theme_times

ggsave(
  filename = file.path(output_plots_dir, "wgcna_soft_threshold_scale_free_fit.png"),
  plot = p_sft_fit,
  width = 8.5,
  height = 6,
  dpi = 300
)

p_sft_conn <- ggplot(fit_df, aes(x = Power, y = MeanConnectivity)) +
  geom_point(size = 3, color = status_colors["IGT"]) +
  geom_line(linewidth = 0.7, color = status_colors["IGT"]) +
  geom_vline(xintercept = soft_power, linetype = "dashed", color = status_colors["T2D"]) +
  geom_text(aes(label = Power), nudge_y = max(fit_df$MeanConnectivity, na.rm = TRUE) * 0.02, size = 3.2, family = "serif") +
  labs(
    title = "Mean connectivity across candidate powers",
    subtitle = paste0("Chosen soft-thresholding power = ", soft_power),
    x = "Soft-thresholding power",
    y = "Mean connectivity"
  ) +
  base_theme_times

ggsave(
  filename = file.path(output_plots_dir, "wgcna_soft_threshold_mean_connectivity.png"),
  plot = p_sft_conn,
  width = 8.5,
  height = 6,
  dpi = 300
)

# Network construction and module detection

cat("Constructing network and detecting modules...\n")

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  corType = "bicor",
  maxPOutliers = 0.1,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  deepSplit = 2,
  verbose = 5
)

module_colors <- net$colors
module_eigengenes <- orderMEs(net$MEs)

cat("Modules detected (including grey): ", length(unique(module_colors)), "\n\n", sep = "")

# Gene annotation and module assignment table

cat("Preparing module assignment table...\n")

gene_ids <- colnames(datExpr)

module_assignment <- tibble(
  gene_id = gene_ids,
  module_color = module_colors
)

if (file.exists(input_gene_annot) && file.info(input_gene_annot)$size > 0) {
  gene_annot <- readr::read_tsv(input_gene_annot, show_col_types = FALSE)

  # gene_annotation.tsv uses original_gene_id as the primary key matching
  # the Ensembl IDs in the expression matrix; join on whichever column matches
  candidate_gene_id_cols <- c("original_gene_id", "ensembl_gene_id", "gene_id")
  annot_gene_col <- candidate_gene_id_cols[candidate_gene_id_cols %in% colnames(gene_annot)][1]

  if (!is.na(annot_gene_col)) {
    gene_annot_small <- gene_annot %>%
      select(all_of(annot_gene_col), any_of(c("gene_symbol", "gene_name", "biotype", "ensembl_gene_id"))) %>%
      distinct(across(all_of(annot_gene_col)), .keep_all = TRUE)

    module_assignment <- module_assignment %>%
      left_join(gene_annot_small, by = setNames(annot_gene_col, "gene_id"))
  }
}

module_counts <- module_assignment %>%
  count(module_color, sort = TRUE, name = "n_genes")

readr::write_csv(module_assignment, file.path(output_modules_dir, "wgcna_gene_module_assignments.csv"))
readr::write_csv(module_counts, file.path(output_modules_dir, "wgcna_module_sizes.csv"))

# Dendrogram with module colours

cat("Creating gene dendrogram and module colour plot...\n")

png(
  filename = file.path(output_modules_dir, "wgcna_gene_dendrogram_with_module_colors.png"),
  width = 4200,
  height = 2400,
  res = 300
)
par(mfrow = c(1, 1), family = "serif", mar = c(6, 5, 4, 2))
plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = module_colors[net$blockGenes[[1]]],
  groupLabels = "Module colour",
  main = "Gene dendrogram and module colours",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

# Module eigengenes and trait matrix

cat("Preparing trait matrix for module-trait analysis...\n")

datTraits_for_cor <- datTraits %>%
  mutate(
    status = factor(as.character(status), levels = c("ND", "IGT", "T2D")),
    status_ND = ifelse(status == "ND", 1, 0),
    status_IGT = ifelse(status == "IGT", 1, 0),
    status_T2D = ifelse(status == "T2D", 1, 0),
    age = make_numeric_trait(age),
    bmi = make_numeric_trait(bmi),
    hba1c = make_numeric_trait(hba1c),
    glycaemia_2h = make_numeric_trait(glycaemia_2h),
    size_factor = make_numeric_trait(size_factor)
  )

trait_cols_for_cor <- c("status_ND", "status_IGT", "status_T2D", "age", "bmi", "hba1c", "glycaemia_2h", "size_factor")

trait_data_numeric <- datTraits_for_cor %>%
  select(all_of(trait_cols_for_cor))

valid_trait_cols <- names(trait_data_numeric)[colSums(!is.na(trait_data_numeric)) >= 10]
trait_data_numeric <- trait_data_numeric %>% select(all_of(valid_trait_cols))

ME_mat <- as.data.frame(module_eigengenes)

module_names <- colnames(ME_mat)
trait_names <- colnames(trait_data_numeric)

cor_mat <- matrix(
  NA_real_,
  nrow = length(module_names),
  ncol = length(trait_names),
  dimnames = list(module_names, trait_names)
)

p_mat <- matrix(
  NA_real_,
  nrow = length(module_names),
  ncol = length(trait_names),
  dimnames = list(module_names, trait_names)
)

for (i in seq_along(module_names)) {
  for (j in seq_along(trait_names)) {
    x <- ME_mat[[module_names[i]]]
    y <- trait_data_numeric[[trait_names[j]]]
    keep <- complete.cases(x, y)

    if (sum(keep) >= 10) {
      tmp <- safe_cor_p(x[keep], y[keep])
      cor_mat[i, j] <- tmp[1]
      p_mat[i, j] <- tmp[2]
    }
  }
}

text_mat <- matrix(
  paste0(sprintf("%.2f", cor_mat), "\n", "p=", formatC(p_mat, format = "e", digits = 2)),
  nrow = nrow(cor_mat),
  dimnames = dimnames(cor_mat)
)

readr::write_csv(
  tibble::as_tibble(cor_mat, rownames = "module"),
  file.path(output_trait_dir, "wgcna_module_trait_correlations.csv")
)

readr::write_csv(
  tibble::as_tibble(p_mat, rownames = "module"),
  file.path(output_trait_dir, "wgcna_module_trait_pvalues.csv")
)

ME_export <- ME_mat %>%
  tibble::rownames_to_column("sampleid")

readr::write_csv(
  ME_export,
  file.path(output_trait_dir, "wgcna_module_eigengenes_by_sample.csv")
)

# Module-trait heatmap

cat("Creating module-trait heatmap...\n")

ht_col_fun <- circlize::colorRamp2(
  c(-1, 0, 1),
  c(status_colors["ND"], "white", status_colors["T2D"])
)

png(
  filename = file.path(output_trait_dir, "wgcna_module_trait_heatmap.png"),
  width = 3200,
  height = 2600,
  res = 300
)

ht <- Heatmap(
  cor_mat,
  name = "Correlation",
  col = ht_col_fun,
  rect_gp = grid::gpar(col = "grey85"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = grid::gpar(fontsize = 10, fontfamily = "Times"),
  column_names_gp = grid::gpar(fontsize = 10, fontfamily = "Times"),
  column_names_rot = 45,
  row_title = "Modules",
  column_title = "Clinical traits",
  row_title_gp = grid::gpar(fontsize = 13, fontfamily = "Times", fontface = "bold"),
  column_title_gp = grid::gpar(fontsize = 13, fontfamily = "Times", fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      text_mat[i, j],
      x, y,
      gp = grid::gpar(fontsize = 8.5, fontfamily = "Times")
    )
  },
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 11, fontfamily = "Times", fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10, fontfamily = "Times")
  )
)

draw(ht)
dev.off()

# Module size bar plot

cat("Creating module size bar plot...\n")

module_counts_plot_df <- module_counts %>%
  mutate(module_color = factor(module_color, levels = module_color))

p_module_sizes <- ggplot(module_counts_plot_df, aes(x = module_color, y = n_genes, fill = module_color)) +
  geom_col(width = 0.75, color = neutral_dark) +
  scale_fill_identity() +
  labs(
    title = "WGCNA module sizes",
    subtitle = "Number of genes assigned to each module",
    x = "Module colour",
    y = "Number of genes"
  ) +
  base_theme_times +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(output_modules_dir, "wgcna_module_sizes_barplot.png"),
  plot = p_module_sizes,
  width = 9,
  height = 6,
  dpi = 300
)

# Module significance summary for disease status

cat("Summarising module significance against disease status...\n")

me_long_summary <- list()

for (me in module_names) {
  for (trait in intersect(c("status_ND", "status_IGT", "status_T2D"), colnames(trait_data_numeric))) {
    keep <- complete.cases(ME_mat[[me]], trait_data_numeric[[trait]])
    if (sum(keep) >= 10) {
      tmp <- safe_cor_p(ME_mat[[me]][keep], trait_data_numeric[[trait]][keep])
      me_long_summary[[length(me_long_summary) + 1]] <- tibble(
        module = me,
        trait = trait,
        correlation = tmp[1],
        p_value = tmp[2],
        significance = p_to_stars(tmp[2])
      )
    }
  }
}

module_trait_summary <- bind_rows(me_long_summary) %>%
  arrange(p_value)

readr::write_csv(
  module_trait_summary,
  file.path(output_trait_dir, "wgcna_module_status_associations_summary.csv")
)

# Save RData objects

cat("Saving RData objects...\n")

save(
  datExpr,
  datTraits,
  gsg,
  sft,
  fit_df,
  soft_power,
  net,
  module_colors,
  module_eigengenes,
  module_assignment,
  module_counts,
  cor_mat,
  p_mat,
  text_mat,
  module_trait_summary,
  file = file.path(output_rdata_dir, "05B_wgcna_network_construction.RData")
)

# summary

cat("\n========================================\n")
cat("05B COMPLETE\n")
cat("========================================\n")
cat("Samples used: ", nrow(datExpr), "\n", sep = "")
cat("Genes used: ", ncol(datExpr), "\n", sep = "")
cat("Chosen soft-thresholding power: ", soft_power, "\n", sep = "")
cat("Modules detected (including grey): ", length(unique(module_colors)), "\n", sep = "")
cat("Outputs written to: ", output_dir, "\n", sep = "")
cat("========================================\n")
