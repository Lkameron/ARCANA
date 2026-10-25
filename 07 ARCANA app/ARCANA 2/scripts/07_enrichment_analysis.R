# 07_enrichment_analysis.R (ARCANA — generalised)
# Reads all settings from config/project.yaml via config_utils.R

source(file.path(Sys.getenv("ARCANA_SCRIPTS_DIR", "scripts"), "config_utils.R"))

# 07_enrichment_analysis.R
#
# Purpose:
#   Functional enrichment analysis (GO Biological Process, Molecular Function,
#   and KEGG pathways) for:
#     (a) DESeq2 differentially expressed gene sets (T2D vs ND, IGT vs ND,
#         T2D vs IGT — both all-gene and lncRNA-only subsets)
#     (b) Genes belonging to each WGCNA focus module
#
# Method:
#   clusterProfiler::enrichGO / enrichKEGG with Ensembl → Entrez conversion
#   via org.Hs.eg.db.  Universe = all tested genes (from DE results).
# Package setup

required_cran <- c(
  "dplyr", "readr", "tibble", "stringr", "ggplot2", "tidyr",
  "forcats", "writexl", "data.table", "scales"
)

required_bioc <- c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
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
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
  # Load the correct OrgDb for the configured species
  orgdb_pkg <- arcana_cfg$orgdb
  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
    BiocManager::install(orgdb_pkg, ask = FALSE, update = FALSE)
  }
  library(orgdb_pkg, character.only = TRUE)
  ORGDB <- get(orgdb_pkg)
})

cat("========================================\n")
cat("07 - Enrichment analysis\n")
cat("========================================\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n", sep = "")

# Paths

project_dir <- getwd()

de_dir  <- file.path(project_dir, "results", "de", "tables")
wgcna_dir <- file.path(project_dir, "results", "wgcna")
annot_path <- file.path(project_dir, "data", "metadata", "gene_annotation.tsv")

out_dir        <- file.path(project_dir, "results", "enrichment")
out_tables_dir <- file.path(out_dir, "tables")
out_plots_dir  <- file.path(out_dir, "plots")
out_rdata_dir  <- file.path(out_dir, "rdata")

for (d in c(out_tables_dir, out_plots_dir, out_rdata_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Plot theme

base_theme <- theme_bw(base_size = 12, base_family = "serif") +
  theme(
    text           = element_text(family = "serif"),
    plot.title     = element_text(face = "bold", size = 13),
    plot.subtitle  = element_text(size = 10),
    axis.title     = element_text(size = 11),
    axis.text      = element_text(size = 9),
    legend.title   = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

status_colors <- c("ND" = "#012345", "IGT" = "#345678", "T2D" = "#8B0000")

# Helper functions

strip_version <- function(x) sub("\\.\\d+$", "", x)

ensembl_to_entrez <- function(ensembl_ids) {
  ids <- strip_version(ensembl_ids)
  mapping <- suppressMessages(
    AnnotationDbi::select(
      ORGDB,
      keys    = ids,
      columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
      keytype = "ENSEMBL"
    )
  )
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  mapping
}

run_go_enrich <- function(gene_entrez, universe_entrez, label,
                           ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2) {
  if (length(gene_entrez) < 5) {
    cat("  Skipping GO (", ont, ") for", label, "— too few genes\n")
    return(NULL)
  }
  tryCatch(
    enrichGO(
      gene          = gene_entrez,
      universe      = universe_entrez,
      OrgDb         = ORGDB,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    ),
    error = function(e) { cat("  GO", ont, "error:", conditionMessage(e), "\n"); NULL }
  )
}

run_kegg_enrich <- function(gene_entrez, universe_entrez, label,
                             p_cutoff = 0.05, q_cutoff = 0.2) {
  if (length(gene_entrez) < 5) {
    cat("  Skipping KEGG for", label, "— too few genes\n")
    return(NULL)
  }
  tryCatch(
    enrichKEGG(
      gene          = gene_entrez,
      universe      = universe_entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    ),
    error = function(e) { cat("  KEGG error:", conditionMessage(e), "\n"); NULL }
  )
}

save_enrich_csv <- function(result_obj, path) {
  if (is.null(result_obj)) return(invisible(NULL))
  df <- as.data.frame(result_obj)
  if (nrow(df) == 0) return(invisible(NULL))
  readr::write_csv(df, path)
  cat("  Saved:", path, "\n")
}

dot_plot_custom <- function(result_obj, title, top_n = 20, fill_color = "#2B6CB0") {
  if (is.null(result_obj)) return(NULL)
  df <- as.data.frame(result_obj)
  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description = forcats::fct_reorder(str_wrap(Description, 40), -log10(p.adjust)),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        if (length(parts) == 2) parts[1] / parts[2] else NA_real_
      })
    )

  ggplot(df, aes(x = GeneRatio_num, y = Description,
                 colour = p.adjust, size = Count)) +
    geom_point() +
    scale_colour_gradient(low = fill_color, high = "#E53E3E",
                          name = "Adj. p-value",
                          labels = scales::scientific) +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(title = title, x = "Gene ratio", y = NULL) +
    base_theme +
    theme(axis.text.y = element_text(size = 8))
}

save_plot <- function(p, path, width = 10, height = 8) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("  Plot saved:", path, "\n")
}

# Load gene annotation 

gene_annot <- NULL
if (file.exists(annot_path) && file.info(annot_path)$size > 0) {
  gene_annot <- readr::read_tsv(annot_path, show_col_types = FALSE)
  cat("Gene annotation loaded:", nrow(gene_annot), "genes\n")
}

# Load DE results and build universe

cat("\nLoading differential expression results...\n")

de_all_files <- setNames(
  lapply(arcana_cfg$comparison_names, function(name) file.path(de_dir, paste0("DESeq2_", name, "_all_genes.tsv"))),
  arcana_cfg$comparison_names
)

de_sig_files <- setNames(
  lapply(arcana_cfg$comparison_names, function(name) file.path(de_dir, paste0("DESeq2_", name, "_significant.tsv"))),
  arcana_cfg$comparison_names
)

de_all  <- lapply(de_all_files,  function(p) if (file.exists(p)) readr::read_tsv(p, show_col_types = FALSE) else NULL)
de_sig  <- lapply(de_sig_files,  function(p) if (file.exists(p)) readr::read_tsv(p, show_col_types = FALSE) else NULL)

# Build universe from the union of all tested genes in any comparison
universe_ensembl <- unique(unlist(lapply(de_all, function(d) if (!is.null(d)) d$gene_id else NULL)))
universe_ensembl <- strip_version(na.omit(universe_ensembl))

cat("Universe size (all tested genes):", length(universe_ensembl), "\n")

universe_map   <- ensembl_to_entrez(universe_ensembl)
universe_entrez <- unique(na.omit(universe_map$ENTREZID))
cat("Universe after Ensembl→Entrez mapping:", length(universe_entrez), "unique Entrez IDs\n\n")

#  DE enrichment — one set per comparison

enrichment_results_de <- list()

for (cmp in names(de_sig)) {
  cat("--- Enrichment for DE genes:", cmp, "---\n")

  sig_df <- de_sig[[cmp]]
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    cat("  No significant DE genes found, skipping.\n")
    next
  }

  sig_ensembl <- strip_version(na.omit(sig_df$gene_id))
  sig_map     <- ensembl_to_entrez(sig_ensembl)
  sig_entrez  <- unique(na.omit(sig_map$ENTREZID))

  cat("  Significant genes:", nrow(sig_df), " | mapped Entrez:", length(sig_entrez), "\n")

  go_bp  <- run_go_enrich(sig_entrez, universe_entrez, cmp, ont = "BP")
  go_mf  <- run_go_enrich(sig_entrez, universe_entrez, cmp, ont = "MF")
  kegg   <- run_kegg_enrich(sig_entrez, universe_entrez, cmp)

  enrichment_results_de[[cmp]] <- list(GO_BP = go_bp, GO_MF = go_mf, KEGG = kegg)

  save_enrich_csv(go_bp,  file.path(out_tables_dir, paste0("GO_BP_", cmp, ".csv")))
  save_enrich_csv(go_mf,  file.path(out_tables_dir, paste0("GO_MF_", cmp, ".csv")))
  save_enrich_csv(kegg,   file.path(out_tables_dir, paste0("KEGG_",  cmp, ".csv")))

  cmp_color <- switch(cmp,
    T2D_vs_ND  = "#8B0000",
    IGT_vs_ND  = "#345678",
    T2D_vs_IGT = "#6B4C9A",
    "#2B6CB0"
  )

  p_go_bp <- dot_plot_custom(go_bp,  paste0("GO BP — ", gsub("_", " ", cmp)), fill_color = cmp_color)
  p_go_mf <- dot_plot_custom(go_mf,  paste0("GO MF — ", gsub("_", " ", cmp)), fill_color = cmp_color)
  p_kegg  <- dot_plot_custom(kegg,   paste0("KEGG — ",  gsub("_", " ", cmp)), fill_color = cmp_color)

  save_plot(p_go_bp, file.path(out_plots_dir, paste0("GO_BP_", cmp, ".png")))
  save_plot(p_go_mf, file.path(out_plots_dir, paste0("GO_MF_", cmp, ".png")))
  save_plot(p_kegg,  file.path(out_plots_dir, paste0("KEGG_",  cmp, ".png")))
}

############################
# SECTION 8: Load WGCNA module assignments
############################

cat("\nLoading WGCNA module data...\n")

module_assign_path <- file.path(wgcna_dir, "modules", "wgcna_gene_module_assignments.csv")
focus_module_path  <- file.path(wgcna_dir, "hubgenes", "tables", "focus_modules_selected_for_hub_analysis.csv")

module_assignments <- readr::read_csv(module_assign_path, show_col_types = FALSE)
focus_modules_tbl  <- readr::read_csv(focus_module_path,  show_col_types = FALSE)
focus_modules      <- focus_modules_tbl$module_color

cat("Focus modules:", paste(focus_modules, collapse = ", "), "\n\n")

############################
# SECTION 9: WGCNA module enrichment
############################

enrichment_results_wgcna <- list()

for (mod in focus_modules) {
  cat("--- Enrichment for WGCNA module:", mod, "---\n")

  mod_genes <- module_assignments %>%
    filter(module_color == mod) %>%
    pull(gene_id)

  mod_ensembl <- strip_version(na.omit(mod_genes))
  mod_map     <- ensembl_to_entrez(mod_ensembl)
  mod_entrez  <- unique(na.omit(mod_map$ENTREZID))

  cat("  Module genes:", length(mod_genes), " | mapped Entrez:", length(mod_entrez), "\n")

  go_bp <- run_go_enrich(mod_entrez, universe_entrez, mod, ont = "BP")
  go_mf <- run_go_enrich(mod_entrez, universe_entrez, mod, ont = "MF")
  kegg  <- run_kegg_enrich(mod_entrez, universe_entrez, mod)

  enrichment_results_wgcna[[mod]] <- list(GO_BP = go_bp, GO_MF = go_mf, KEGG = kegg)

  save_enrich_csv(go_bp, file.path(out_tables_dir, paste0("GO_BP_module_", mod, ".csv")))
  save_enrich_csv(go_mf, file.path(out_tables_dir, paste0("GO_MF_module_", mod, ".csv")))
  save_enrich_csv(kegg,  file.path(out_tables_dir, paste0("KEGG_module_",  mod, ".csv")))

  p_go_bp <- dot_plot_custom(go_bp,  paste0("GO BP — ", mod, " module"), fill_color = mod)
  p_go_mf <- dot_plot_custom(go_mf,  paste0("GO MF — ", mod, " module"), fill_color = mod)
  p_kegg  <- dot_plot_custom(kegg,   paste0("KEGG — ",  mod, " module"), fill_color = mod)

  save_plot(p_go_bp, file.path(out_plots_dir, paste0("GO_BP_module_", mod, ".png")))
  save_plot(p_go_mf, file.path(out_plots_dir, paste0("GO_MF_module_", mod, ".png")))
  save_plot(p_kegg,  file.path(out_plots_dir, paste0("KEGG_module_",  mod, ".png")))
}

############################
# SECTION 10: Comparative dotplot across focus modules (GO BP)
############################

cat("\nBuilding cross-module GO summary plot...\n")

top5_per_module <- lapply(names(enrichment_results_wgcna), function(mod) {
  go_bp <- enrichment_results_wgcna[[mod]]$GO_BP
  if (is.null(go_bp)) return(NULL)
  df <- as.data.frame(go_bp)
  if (nrow(df) == 0) return(NULL)
  df %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    mutate(module = mod)
}) %>% bind_rows()

if (nrow(top5_per_module) > 0) {
  top5_per_module <- top5_per_module %>%
    mutate(
      label = str_wrap(Description, 35),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        if (length(parts) == 2) parts[1] / parts[2] else NA_real_
      })
    )

  p_cross <- ggplot(top5_per_module,
                    aes(x = module, y = reorder(label, -log10(p.adjust)),
                        size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_gradient(low = "grey70", high = "#8B0000",
                          name = "-log10\nadj. p-value") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(
      title = "Top GO Biological Process terms per WGCNA focus module",
      subtitle = "Top 5 terms per module, ranked by adjusted p-value",
      x = "Module colour",
      y = NULL
    ) +
    base_theme +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(angle = 30, hjust = 1))

  save_plot(p_cross,
            file.path(out_plots_dir, "GO_BP_cross_module_summary.png"),
            width = 13, height = max(6, nrow(top5_per_module) * 0.35 + 3))
}

############################
# SECTION 11: DE comparison summary enrichment bar chart
############################

cat("Building DE enrichment comparison summary...\n")

de_top_terms <- lapply(names(enrichment_results_de), function(cmp) {
  go_bp <- enrichment_results_de[[cmp]]$GO_BP
  if (is.null(go_bp)) return(NULL)
  df <- as.data.frame(go_bp)
  if (nrow(df) == 0) return(NULL)
  df %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    mutate(comparison = gsub("_", " ", cmp))
}) %>% bind_rows()

if (nrow(de_top_terms) > 0) {
  de_top_terms <- de_top_terms %>%
    mutate(
      label  = str_wrap(Description, 35),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        if (length(parts) == 2) parts[1] / parts[2] else NA_real_
      }),
      comp_color = case_when(
        grepl("T2D vs ND", comparison)  ~ "#8B0000",
        grepl("IGT vs ND", comparison)  ~ "#345678",
        grepl("T2D vs IGT", comparison) ~ "#6B4C9A",
        TRUE ~ "grey50"
      )
    )

  p_de_summary <- ggplot(de_top_terms,
                          aes(x = GeneRatio_num, y = reorder(label, GeneRatio_num),
                              fill = comparison, size = Count)) +
    geom_point(shape = 21, colour = "grey30") +
    scale_fill_manual(values = c(
      "T2D vs ND"  = "#8B0000",
      "IGT vs ND"  = "#345678",
      "T2D vs IGT" = "#6B4C9A"
    ), name = "Comparison") +
    scale_size_continuous(name = "Gene count", range = c(3, 9)) +
    labs(
      title    = "Top GO Biological Process terms across DE comparisons",
      subtitle = "Top 5 terms per comparison by adjusted p-value",
      x = "Gene ratio",
      y = NULL
    ) +
    base_theme +
    theme(axis.text.y = element_text(size = 8))

  save_plot(p_de_summary,
            file.path(out_plots_dir, "GO_BP_DE_comparisons_summary.png"),
            width = 12, height = 8)
}

############################
# SECTION 12: Compile Excel summary
############################

cat("\nCompiling Excel summary...\n")

all_sheets <- list()

for (cmp in names(enrichment_results_de)) {
  go_bp_df <- if (!is.null(enrichment_results_de[[cmp]]$GO_BP)) {
    as.data.frame(enrichment_results_de[[cmp]]$GO_BP)
  } else data.frame()
  kegg_df <- if (!is.null(enrichment_results_de[[cmp]]$KEGG)) {
    as.data.frame(enrichment_results_de[[cmp]]$KEGG)
  } else data.frame()

  sheet_name_go   <- paste0("DE_GOBP_", substr(cmp, 1, 20))
  sheet_name_kegg <- paste0("DE_KEGG_", substr(cmp, 1, 20))

  if (nrow(go_bp_df) > 0) all_sheets[[sheet_name_go]]   <- go_bp_df
  if (nrow(kegg_df)  > 0) all_sheets[[sheet_name_kegg]] <- kegg_df
}

for (mod in names(enrichment_results_wgcna)) {
  go_bp_df <- if (!is.null(enrichment_results_wgcna[[mod]]$GO_BP)) {
    as.data.frame(enrichment_results_wgcna[[mod]]$GO_BP)
  } else data.frame()
  kegg_df <- if (!is.null(enrichment_results_wgcna[[mod]]$KEGG)) {
    as.data.frame(enrichment_results_wgcna[[mod]]$KEGG)
  } else data.frame()

  sheet_name_go   <- paste0("Mod_GOBP_", substr(mod, 1, 20))
  sheet_name_kegg <- paste0("Mod_KEGG_", substr(mod, 1, 20))

  if (nrow(go_bp_df) > 0) all_sheets[[sheet_name_go]]   <- go_bp_df
  if (nrow(kegg_df)  > 0) all_sheets[[sheet_name_kegg]] <- kegg_df
}

if (length(all_sheets) > 0) {
  writexl::write_xlsx(all_sheets,
    file.path(project_dir, "GSE164416_enrichment_analysis_summary.xlsx"))
  cat("Excel summary written.\n")
}

############################
# SECTION 13: Save RData
############################

save(
  enrichment_results_de,
  enrichment_results_wgcna,
  universe_entrez,
  universe_map,
  focus_modules,
  file = file.path(out_rdata_dir, "07_enrichment.RData")
)

############################
# SECTION 14: Final report
############################

cat("\n========================================\n")
cat("07 COMPLETE\n")
cat("========================================\n")
cat("DE comparisons enriched:", length(enrichment_results_de), "\n")
cat("WGCNA modules enriched:", length(enrichment_results_wgcna), "\n")
cat("Outputs written to: results/enrichment/\n")
cat("  tables/ — CSV enrichment result tables\n")
cat("  plots/  — dot plots per comparison and per module\n")
cat("  rdata/07_enrichment.RData\n")
cat("  GSE164416_enrichment_analysis_summary.xlsx\n")
cat("========================================\n")
