# SCRIPT: 05A_prepare_wgcna_inputs.R
#
# WHAT THIS SCRIPT DOES:
#   1. Loads the VST-ready expression matrix for WGCNA
#   2. Loads metadata and clinical traits
#   3. Matches sample order across all inputs
#   4. Runs basic WGCNA input QC with goodSamplesGenes()
#   5. Produces sample clustering plots for outlier inspection
#   6. Writes cleaned expression and trait files for downstream WGCNA scripts

# Package setup

required_cran <- c(
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tibble",
  "ggplot2"
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
library(stringr)
library(tibble)
library(ggplot2)
library(WGCNA)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

cat("========================================\n")
cat("05A - Prepare WGCNA inputs\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n\n", sep = "")

# Define paths

project_dir <- getwd()

input_expr_file <- file.path(project_dir, "data", "processed", "qc", "normalized_counts_vst_input.tsv")
input_traits_file <- file.path(project_dir, "data", "metadata", "clinical_traits.csv")
fallback_traits_file <- file.path(project_dir, "data", "processed", "qc", "sample_metadata_ND_IGT_T2D.csv")

output_dir <- file.path(project_dir, "results", "wgcna")
output_tables <- file.path(output_dir, "tables")
output_plots <- file.path(output_dir, "plots")
output_rdata <- file.path(output_dir, "rdata")

dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(output_rdata, recursive = TRUE, showWarnings = FALSE)


check_file_exists_nonempty <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) {
    stop("File exists but is empty: ", path)
  }
}

standardise_ids <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x
}

find_sample_column <- function(df, df_name) {
  possible_names <- c(
    "sampleid", "sample_id", "sample", "SampleID", "Sample",
    "geo_accession", "GEO_accession", "geoaccession", "donor_id", "donorid"
  )
  hit <- possible_names[possible_names %in% colnames(df)]
  if (length(hit) == 0) {
    stop(
      "Could not find a sample ID column in ", df_name, ".\n",
      "Please add one of these column names: ",
      paste(possible_names, collapse = ", ")
    )
  }
  hit[1]
}

extract_optional_numeric <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[hit[1]]]))
}

extract_optional_character <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
  as.character(df[[hit[1]]])
}

# Check files

check_file_exists_nonempty(input_expr_file)

if (file.exists(input_traits_file) && file.info(input_traits_file)$size > 0) {
  traits_source <- input_traits_file
} else if (file.exists(fallback_traits_file) && file.info(fallback_traits_file)$size > 0) {
  traits_source <- fallback_traits_file
  cat("clinical_traits.csv missing or empty; using fallback processed QC metadata file instead.\n\n")
} else {
  stop(
    "No usable trait file found. Expected one of:\n",
    input_traits_file, "\n",
    fallback_traits_file
  )
}

cat("All required input files found.\n\n")

# Read input data

cat("Reading expression and trait files...\n")

expr_df <- data.table::fread(input_expr_file, data.table = FALSE)
traits_raw <- readr::read_csv(traits_source, show_col_types = FALSE)

cat("Expression matrix dimensions (raw file): ", nrow(expr_df), " rows x ", ncol(expr_df), " columns\n", sep = "")
cat("Clinical traits dimensions: ", nrow(traits_raw), " rows x ", ncol(traits_raw), " columns\n", sep = "")
cat("Traits source: ", traits_source, "\n\n", sep = "")

# Prepare expression matrix

cat("Preparing expression matrix for WGCNA...\n")

gene_id_col <- colnames(expr_df)[1]

expr_mat <- expr_df %>%
  tibble::column_to_rownames(var = gene_id_col) %>%
  as.matrix()

mode(expr_mat) <- "numeric"

datExpr0 <- as.data.frame(t(expr_mat))
rownames(datExpr0) <- standardise_ids(rownames(datExpr0))

cat("WGCNA expression dimensions after transpose: ",
    nrow(datExpr0), " samples x ", ncol(datExpr0), " genes\n\n", sep = "")

# Prepare trait table

traits_sample_col <- find_sample_column(traits_raw, basename(traits_source))
cat("Traits sample column detected: ", traits_sample_col, "\n\n", sep = "")

traits_df <- tibble(
  sampleid = standardise_ids(traits_raw[[traits_sample_col]]),
  status = extract_optional_character(traits_raw, c("status")),
  age = extract_optional_numeric(traits_raw, c("age")),
  sex = extract_optional_character(traits_raw, c("sex")),
  bmi = extract_optional_numeric(traits_raw, c("bmi")),
  hba1c = extract_optional_numeric(traits_raw, c("hba1c")),
  glycaemia_2h = extract_optional_numeric(traits_raw, c("glycaemia_2h", "glycaemia2h")),
  size_factor = extract_optional_numeric(traits_raw, c("size_factor", "sizefactor"))
) %>%
  distinct(sampleid, .keep_all = TRUE) %>%
  filter(!is.na(sampleid), sampleid != "")

traits_df <- traits_df %>%
  mutate(
    status = case_when(
      status %in% c("IFG", "Impaired Fasting Glucose") ~ "IGT",
      TRUE ~ status
    ),
    status = factor(status, levels = c("ND", "IGT", "T2D"))
  )

# Match samples

cat("Matching expression samples to trait table...\n")

expr_samples <- rownames(datExpr0)
trait_samples <- traits_df$sampleid

samples_in_expr_not_traits <- setdiff(expr_samples, trait_samples)
samples_in_traits_not_expr <- setdiff(trait_samples, expr_samples)
shared_samples <- intersect(expr_samples, trait_samples)

cat("Samples in expression only: ", length(samples_in_expr_not_traits), "\n", sep = "")
cat("Samples in traits only: ", length(samples_in_traits_not_expr), "\n", sep = "")
cat("Shared samples: ", length(shared_samples), "\n\n", sep = "")

if (length(shared_samples) < 10) {
  stop("Too few overlapping samples between expression and trait files. Please check sample naming.")
}

datExpr1 <- datExpr0[shared_samples, , drop = FALSE]
datTraits1 <- traits_df %>%
  filter(sampleid %in% shared_samples) %>%
  dplyr::slice(match(shared_samples, sampleid))

stopifnot(all(rownames(datExpr1) == datTraits1$sampleid))

# WGCNA QC

cat("Running goodSamplesGenes() QC...\n")

gsg <- WGCNA::goodSamplesGenes(datExpr1, verbose = 3)

removed_genes_df <- tibble(
  gene_id = colnames(datExpr1),
  removed_by_goodSamplesGenes = !gsg$goodGenes
)

removed_samples_df <- tibble(
  sampleid = rownames(datExpr1),
  removed_by_goodSamplesGenes = !gsg$goodSamples,
  reason = ifelse(!gsg$goodSamples, "Failed goodSamplesGenes", "")
)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing ", sum(!gsg$goodGenes), " genes flagged by goodSamplesGenes.\n", sep = "")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing ", sum(!gsg$goodSamples), " samples flagged by goodSamplesGenes.\n", sep = "")
  }

  datExpr1 <- datExpr1[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  datTraits1 <- datTraits1 %>%
    filter(sampleid %in% rownames(datExpr1)) %>%
    dplyr::slice(match(rownames(datExpr1), sampleid))
} else {
  cat("All samples and genes passed goodSamplesGenes().\n")
}

cat("\nRetained after QC: ", nrow(datExpr1), " samples x ", ncol(datExpr1), " genes\n\n", sep = "")

# Sample clustering plot 

cat("Creating sample clustering plot...\n")

status_colors <- c(
  "ND"  = "#012345",  # Midnight Navy
  "IGT" = "#345678",  # Deep Sky
  "T2D" = "#8b0000"   # Dark Crimson
)

# Build a label vector 
if ("donor_id" %in% colnames(datTraits1)) {
  label_vec <- paste0(datTraits1$donor_id, "_", datTraits1$status)
} else {
  label_vec <- paste0(datTraits1$sampleid, "_", datTraits1$status)
}


status_simple <- as.character(datTraits1$status)

# Dendrogram object
sampleTree <- hclust(dist(datExpr1), method = "average")

# Plot with larger margins, Times, and coloured labels
png(
  filename = file.path(output_plots, "sample_clustering_before_outlier_removal.png"),
  width = 3600, height = 2200, res = 300
)
par(
  mar = c(10, 5, 5, 2),        
  family = "Times",
  cex = 0.7                    
)

plot(
  sampleTree,
  labels = label_vec,
  main  = "Sample clustering before outlier removal",
  sub   = "",
  xlab  = "",
  ylab  = "Height",
  cex.main = 1.4,
  cex.lab  = 1.2,
  cex.axis = 1.0
)

# Colour the labels by status 
# Map label positions back to samples:
label_order <- sampleTree$labels[sampleTree$order]
status_ordered <- status_simple[match(label_order, label_vec)]
label_cols <- status_colors[status_ordered]


axis(
  side = 1,
  at   = seq_along(label_order),
  labels = FALSE
)

text(
  x = seq_along(label_order),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = label_order,
  srt = 90,
  xpd = TRUE,
  adj = 1,
  cex = 0.6,
  col = label_cols
)

dev.off()

# Sample distribution bar plot (improved)

if (any(!is.na(datTraits1$status))) {
  trait_df <- datTraits1 %>%
    mutate(status = factor(as.character(status),
                           levels = c("ND", "IGT", "T2D")))

  status_colors <- c(
    "ND"  = "#012345",
    "IGT" = "#345678",
    "T2D" = "#8b0000"
  )

  # Simple chi‑square test of equal proportions
  tab_counts <- table(trait_df$status)
  chisq_res  <- chisq.test(tab_counts)
  pval       <- chisq_res$p.value

  sig_label <- dplyr::case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ "ns"
  )

  # Y‑position for significance label
  y_max <- max(tab_counts)
  y_sig <- y_max * 1.08

  p_status <- ggplot(trait_df, aes(x = status, fill = status)) +
    geom_bar(width = 0.7, color = "#2B1B17", alpha = 0.9) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    labs(
      title    = "WGCNA input sample distribution",
      subtitle = paste0(
        "Samples retained for co-expression analysis (chi-square p = ",
        signif(pval, 3), ", ", sig_label, ")"
      ),
      x = "Disease status",
      y = "Number of samples"
    ) +
    theme_bw(base_size = 12, base_family = "serif") +
    theme(
      text            = element_text(family = "Times"),
      plot.title      = element_text(face = "bold", size = 14),
      plot.subtitle   = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    # significance label centred above all three bars
    annotate(
      "text",
      x = 2,           
      y = y_sig,
      label = sig_label,
      size = 5,
      family = "Times"
    ) +
    coord_cartesian(ylim = c(0, y_sig * 1.1), clip = "off")

  ggsave(
    filename = file.path(output_plots, "wgcna_sample_distribution_by_status.png"),
    plot     = p_status,
    width    = 7,
    height   = 5.5,
    dpi      = 300
  )
}

# Save 

cat("Writing cleaned WGCNA input files...\n")

wgcna_expression_output <- datExpr1 %>%
  rownames_to_column("sampleid")

readr::write_tsv(
  wgcna_expression_output,
  file.path(output_tables, "wgcna_expression_input.tsv")
)

readr::write_csv(
  datTraits1,
  file.path(output_tables, "wgcna_traits_input.csv")
)

readr::write_tsv(
  removed_genes_df,
  file.path(output_tables, "wgcna_removed_genes.tsv")
)

readr::write_tsv(
  removed_samples_df,
  file.path(output_tables, "wgcna_removed_samples.tsv")
)

save(
  datExpr1,
  datTraits1,
  sampleTree,
  removed_genes_df,
  removed_samples_df,
  file = file.path(output_rdata, "05A_wgcna_inputs.RData")
)

# summary

cat("\n========================================\n")
cat("05A COMPLETE\n")
cat("========================================\n")
cat("Expression file used: ", input_expr_file, "\n", sep = "")
cat("Traits file used: ", traits_source, "\n\n", sep = "")
cat("Final WGCNA matrix: ", nrow(datExpr1), " samples x ", ncol(datExpr1), " genes\n", sep = "")
cat("Outputs written to: ", output_dir, "\n", sep = "")
cat("========================================\n")