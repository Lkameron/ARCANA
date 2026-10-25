#!/usr/bin/env Rscript
# N1_nested_split.R <RUN_ID> <SEED>
#
# Purpose:
#   Create ONE stratified train/held-out split BEFORE any WGCNA module
#   detection or hub-gene selection touches the held-out samples. This is the
#   leakage fix for the ML classification claim in script 08: hub-gene
#   features were originally selected using module-trait correlation across
#   ALL 98 samples (including the labels of every "held-out" CV fold), which
#   inflates the reported classification AUC. Re-deriving WGCNA + hub genes
#   from TRAIN samples only, then evaluating once on the never-touched
#   held-out set, removes that leakage.
#
#   Run multiple times with different RUN_ID/SEED to check whether the
#   held-out AUC is stable across different random splits, rather than
#   trusting a single split (with only ~20 held-out samples, a single split
#   can easily be lucky or unlucky).
#
# Outputs (per run):
#   results/nested_validation/<RUN_ID>/tables/train_test_split.csv
#   results/nested_validation/<RUN_ID>/tables/wgcna_expression_input_TRAIN.tsv
#   results/nested_validation/<RUN_ID>/tables/wgcna_traits_input_TRAIN.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) args[1] else "run1"
SEED   <- if (length(args) >= 2) as.integer(args[2]) else 42

set.seed(SEED)
cat("RUN_ID =", RUN_ID, " | SEED =", SEED, "\n\n")

project_dir <- getwd()
out_dir <- file.path(project_dir, "results", "nested_validation", RUN_ID, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- read_csv(file.path(project_dir, "results", "wgcna", "tables", "wgcna_traits_input.csv"),
                    show_col_types = FALSE)
expr   <- read_tsv(file.path(project_dir, "results", "wgcna", "tables", "wgcna_expression_input.tsv"),
                    show_col_types = FALSE)

cat("Total samples:", nrow(traits), "\n")
print(table(traits$status))

# Stratified 80/20 split by status
split_list <- lapply(split(traits$sampleid, traits$status), function(ids) {
  n_test <- max(1, round(length(ids) * 0.2))
  test_ids <- sample(ids, n_test)
  data.frame(sampleid = ids, split = ifelse(ids %in% test_ids, "test", "train"))
})
split_df <- bind_rows(split_list)

cat("\nSplit counts:\n")
print(table(traits$status[match(split_df$sampleid, traits$sampleid)], split_df$split))

write_csv(split_df, file.path(out_dir, "train_test_split.csv"))

train_ids <- split_df$sampleid[split_df$split == "train"]
test_ids  <- split_df$sampleid[split_df$split == "test"]

cat("\nTrain n =", length(train_ids), " | Test n =", length(test_ids), "\n")

traits_train <- traits %>% filter(sampleid %in% train_ids)
expr_train   <- expr %>% filter(sampleid %in% train_ids)

write_csv(traits_train, file.path(out_dir, "wgcna_traits_input_TRAIN.csv"))
write_tsv(expr_train, file.path(out_dir, "wgcna_expression_input_TRAIN.tsv"))

cat("\nWritten to", out_dir, "\n")
cat("\nN1 COMPLETE (", RUN_ID, ")\n", sep = "")
