# 00_download_data.R
# the actions performed by this script are:
# 1. Downloading raw count matrix from GEO
# 2. Extracting sample metadata from GEO Series Matrix
# 3. Integrating clinical metadata from paper supplementary files
# 4. Annotating genes with symbols, names, and biotypes from org.Hs.eg.db and GENCODE
# 5. Exporting cleaned datasets for downstream QC and normalization


# Loading the CRAN packages
required_packages <- c(
  "data.table",    
  "dplyr",         
  "readr",         
  "stringr",       
  "tibble",        
  "writexl",       
  "readxl",        
  "R.utils"        
)

# Bioconductor packages
bioc_packages <- c(
  "GEOquery",      
  "AnnotationDbi", 
  "org.Hs.eg.db"   
)

# Installing BiocManager code
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

# libraries
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(writexl)
library(readxl)
library(R.utils)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)

# For reproducibility purposes, print session info
cat("R version:", R.version.string, "\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")


# Set working directory 
project_dir <- getwd()

# Creating subdirectories for better organisation of data files
raw_dir <- file.path(project_dir, "data", "raw", "GSE164416")
metadata_dir <- file.path(project_dir, "data", "metadata")
processed_dir <- file.path(project_dir, "data", "processed")


dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# file paths for downloaded data
count_gz <- file.path(raw_dir, "GSE164416_DP_htseq_counts.txt.gz")
count_txt <- file.path(raw_dir, "GSE164416_DP_htseq_counts.txt")
gtf_gz <- file.path(raw_dir, "gencode.v46.annotation.gtf.gz")

# Path to clinical metadata from paper supplement (Supplementary Table 1)
# Downloaded from: https://www.nature.com/articles/s42255-021-00420-9
clinical_xlsx <- file.path(project_dir, "42255_2021_420_MOESM4_ESM.xlsx")

cat("Project directories created successfully.\n\n")


# GSE164416 supplementary file URL
count_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/",
  "GSE164416/suppl/GSE164416_DP_htseq_counts.txt.gz"
)

# Downloading count matrix t
if (!file.exists(count_gz)) {
  cat("Downloading count matrix from GEO FTP...\n")
  download.file(count_url, count_gz, mode = "wb")
  cat("Download complete.\n")
} else {
  cat("Count matrix already downloaded.\n")
}


if (!file.exists(count_txt)) {
  cat("Decompressing count matrix...\n")
  R.utils::gunzip(count_gz, remove = FALSE, overwrite = TRUE)
  cat("Decompression complete.\n\n")
} else {
  cat("Count matrix already decompressed.\n\n")
}


# GEO Series Matrix file containing sample annotations
cat("Downloading GEO Series Matrix metadata...\n")
gse_list <- GEOquery::getGEO(
  "GSE164416",
  GSEMatrix = TRUE,
  destdir = raw_dir
)

# phenotype 
gse_object <- gse_list[[1]]
geo_pheno <- Biobase::pData(gse_object)

cat("GEO metadata downloaded successfully.\n")
cat("Number of samples in GEO:", nrow(geo_pheno), "\n\n")

# next i need to load and process the count matrix [24/02/2026]
cat("Loading count matrix...\n")

counts_df <- data.table::fread(count_txt)

gene_id_col <- colnames(counts_df)[1]
count_matrix <- counts_df %>%
  tibble::column_to_rownames(var = gene_id_col) %>%
  as.matrix()

mode(count_matrix) <- "numeric"

sample_names <- colnames(count_matrix)
gene_ids <- rownames(count_matrix)

cat("Count matrix loaded successfully.\n")
cat("Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n\n")

# need to find the supplementary file with clinical metadata for the donors. The paper mentions that the clinical characteristics of the donors are provided in Supplementary Table 1, which is available as an Excel file. I will check if this file exists in the project directory and load it.
# found supplementary file but have to pay for a nature subscription to access it. # nolint
# SUPPLEMENTARY FILE DOWNLOADED 

cat("Loading clinical metadata from supplementary file...\n")

if (!file.exists(clinical_xlsx)) {
  stop(paste0(
    "\nERROR: Clinical metadata file not found!\n",
    "Expected location: ", clinical_xlsx, "\n",
    "Please download '42255_2021_420_MOESM4_ESM-3.xlsx' from:\n",
    "https://www.nature.com/articles/s42255-021-00420-9#Sec23\n",
    "And place it in: ", project_dir, "\n"
  ))
}

clinical_data <- readxl::read_excel(
  clinical_xlsx,
  sheet = 1,
  col_names = TRUE
) %>%
  # Rename columns 
  dplyr::rename(
    donor_id = 1,
    diabetes_status_clinical = 2,
    age_clinical = 3,                
    glycaemia_2h_ogtt = 4,           
    hba1c_clinical = 5,              
    bmi_clinical = 6                 
  ) %>%
  dplyr::mutate(
    # Clean donor ID 
    donor_id = stringr::str_trim(donor_id),
    
    # Standardise diabetes labels
    # Note: Supplementary file has "T3D" instead of "T3cD"
    status_clinical = dplyr::case_when(
      diabetes_status_clinical == "ND" ~ "ND",
      diabetes_status_clinical == "IGT" ~ "IGT",
      diabetes_status_clinical == "T2D" ~ "T2D",
      diabetes_status_clinical %in% c("T3D", "T3cD") ~ "T3cD",
      TRUE ~ "Unknown"
    ),
    
    
    
    age = as.numeric(age_clinical),
    bmi = as.numeric(bmi_clinical),
    hba1c = as.numeric(hba1c_clinical),
    glycaemia_2h = as.numeric(glycaemia_2h_ogtt)
  ) %>%
  # selected columns for merging data 
  dplyr::select(
    donor_id,
    status_clinical,
    age,
    bmi,
    hba1c,
    glycaemia_2h
  )

cat("Clinical data loaded for", nrow(clinical_data), "donors.\n\n")


cat("Building integrated sample metadata...\n")

sample_metadata <- tibble(
  sample = sample_names,
  donor_id = stringr::str_extract(sample_names, "DP\\d+")
)

# Extract relevant fields from GEO phenotype data
geo_metadata <- geo_pheno %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    # Extract donor ID from title field
    donor_id = stringr::str_extract(title, "DP\\d+"),
    
    # Parse diabetes status
    diabetes_status_geo = stringr::str_remove(
      characteristics_ch1,
      "^diabetes status:\\s*"
    ),
    
    # Standardise status labels
    status_geo = dplyr::case_when(
      diabetes_status_geo == "ND" ~ "ND",
      diabetes_status_geo == "IGT" ~ "IGT",
      diabetes_status_geo == "T2D" ~ "T2D",
      diabetes_status_geo == "T3cD" ~ "T3cD",
      TRUE ~ "Unknown"
    )
  ) %>%
  
  dplyr::select(
    donor_id,
    title,
    geo_accession,
    status_geo
  )

# Merging all metadata sources by donor ID
# Priority: clinical data > GEO data 
sample_metadata <- sample_metadata %>%
  dplyr::left_join(geo_metadata, by = "donor_id") %>%
  dplyr::left_join(clinical_data, by = "donor_id") %>%
  dplyr::mutate(
    
    status = dplyr::coalesce(status_clinical, status_geo),
    
    # NOTE: Sex is not available in the supplementary file
    sex = NA_character_
  ) %>%
  # Select final column order
  dplyr::select(
    sample,
    donor_id,
    geo_accession,
    title,
    status,
    age,
    sex,
    bmi,
    hba1c,
    glycaemia_2h
  )

# Print summary statistics so that when i run the script it gives all data at the end for easy methods writings 
cat("\n========================================\n")
cat("Sample Metadata Summary\n")
cat("========================================\n")
cat("Total samples:", nrow(sample_metadata), "\n")
cat("\nDiabetes status distribution:\n")
cat("  ND samples:", sum(sample_metadata$status == "ND", na.rm = TRUE), "\n")
cat("  IGT samples:", sum(sample_metadata$status == "IGT", na.rm = TRUE), "\n")
cat("  T2D samples:", sum(sample_metadata$status == "T2D", na.rm = TRUE), "\n")
cat("  T3cD samples:", sum(sample_metadata$status == "T3cD", na.rm = TRUE), "\n")
cat("  Unknown:", sum(sample_metadata$status == "Unknown", na.rm = TRUE), "\n")
cat("\nClinical variables available:\n")
cat("  Age:", sum(!is.na(sample_metadata$age)), "samples\n")
cat("  Sex:", sum(!is.na(sample_metadata$sex)), "samples (not in supplementary data)\n")
cat("  BMI:", sum(!is.na(sample_metadata$bmi)), "samples\n")
cat("  HbA1c:", sum(!is.na(sample_metadata$hba1c)), "samples\n")
cat("========================================\n\n")


# Map Ensembl gene IDs to gene symbols and full names using org.Hs.eg.db as ensemble keeps crashing i have to use the offline version. 
# [04/03/2026]

cat("Annotating genes with symbols and names...\n")

# Remove version numbers from Ensembl IDs (e.g., ENSG00000000003.15 -> ENSG00000000003)
gene_ids_clean <- stringr::str_remove(gene_ids, "\\.[0-9]+$")

# initial annotation table
gene_annotation <- tibble(
  original_gene_id = gene_ids,        
  ensembl_gene_id = gene_ids_clean    
) %>%
  dplyr::mutate(
    
    gene_symbol = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_gene_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"  
    ) %>% unname(),
    
    
    gene_name = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_gene_id,
      column = "GENENAME",
      keytype = "ENSEMBL",
      multiVals = "first"
    ) %>% unname()
  ) %>%
  dplyr::mutate(
    
    gene_symbol = dplyr::if_else(is.na(gene_symbol), ensembl_gene_id, gene_symbol)
  )

cat("Gene symbol mapping complete.\n")
cat("Mapped", sum(!is.na(gene_annotation$gene_symbol)), "genes to symbols.\n\n")


# gene biotype classification 
# used GENCODE v46 GTF for this 

cat("Downloading GENCODE annotation for gene biotypes...\n")

gtf_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
)

if (!file.exists(gtf_gz)) {
  cat("Downloading GTF from GENCODE FTP...\n")
  download.file(gtf_url, gtf_gz, mode = "wb")
  cat("Download complete.\n")
}

cat("Parsing GTF file for gene entries...\n")
cat("(This may take 1-2 minutes for ~60,000 genes)\n\n")

gtf_lines <- readr::read_lines(gtf_gz, progress = FALSE)

gtf_lines <- gtf_lines[!stringr::str_starts(gtf_lines, "#")]

cat("Total GTF lines:", format(length(gtf_lines), big.mark = ","), "\n")

gtf_split <- stringr::str_split(gtf_lines, "\t", n = 9)

is_gene <- sapply(gtf_split, function(x) {
  length(x) >= 9 && x[3] == "gene"
})

gene_lines <- gtf_lines[is_gene]

cat("Found", format(length(gene_lines), big.mark = ","), "gene entries\n")


# Attributes format: gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1";
gtf_gene_df <- tibble(line = gene_lines) %>%
  dplyr::mutate(
    
    ensembl_gene_id_versioned = stringr::str_match(
      line,
      'gene_id "([^"]+)"'
    )[, 2],
    
   
    gene_name_gtf = stringr::str_match(
      line,
      'gene_name "([^"]+)"'
    )[, 2],
    
    
    biotype = stringr::str_match(
      line,
      'gene_type "([^"]+)"'
    )[, 2]
  ) %>%
  dplyr::mutate(
    
    ensembl_gene_id = stringr::str_remove(ensembl_gene_id_versioned, "\\.[0-9]+$")
  ) %>%
  dplyr::filter(!is.na(ensembl_gene_id)) %>%
  dplyr::select(ensembl_gene_id, gene_name_gtf, biotype) %>%
  
  dplyr::distinct()

cat("Parsed", format(nrow(gtf_gene_df), big.mark = ","), "unique genes with biotype annotation\n")

# Merge biotype information into gene annotation table
gene_annotation <- gene_annotation %>%
  dplyr::left_join(gtf_gene_df, by = "ensembl_gene_id") %>%
  dplyr::mutate(
    
    gene_name = dplyr::if_else(
      is.na(gene_name) & !is.na(gene_name_gtf),
      gene_name_gtf,
      gene_name
    )
  ) %>%
  # column selection
  dplyr::select(
    original_gene_id,
    ensembl_gene_id,
    gene_symbol,
    gene_name,
    biotype
  )

# Print summary statistics
cat("\n========================================\n")
cat("Gene Annotation Summary\n")
cat("========================================\n")
cat("Total genes:", format(nrow(gene_annotation), big.mark = ","), "\n")
cat("Genes with symbols:", format(sum(!is.na(gene_annotation$gene_symbol)), big.mark = ","), "\n")
cat("Genes with biotype:", format(sum(!is.na(gene_annotation$biotype)), big.mark = ","), "\n")
cat("\nBiotype distribution:\n")
cat("  Protein-coding:", format(sum(gene_annotation$biotype == "protein_coding", na.rm = TRUE), big.mark = ","), "\n")
cat("  lncRNA:", format(sum(gene_annotation$biotype == "lncRNA", na.rm = TRUE), big.mark = ","), "\n")
cat("  Pseudogenes:", format(sum(stringr::str_detect(gene_annotation$biotype, "pseudogene"), na.rm = TRUE), big.mark = ","), "\n")
cat("========================================\n\n")


cat("Exporting processed data files...\n")

# Export count matrix (TSV)
readr::write_tsv(
  counts_df,
  file.path(processed_dir, "counts_matrix.tsv")
)

# Export sample metadata (CSV)
readr::write_csv(
  sample_metadata,
  file.path(metadata_dir, "sample_metadata.csv")
)

# Export gene annotation (TSV )
readr::write_tsv(
  gene_annotation,
  file.path(metadata_dir, "gene_annotation.tsv")
)


summary_table <- tibble(
  metric = c(
    "total_genes",
    "total_samples",
    "nd_samples",
    "igt_samples",
    "t2d_samples",
    "t3cd_samples",
    "age_available",
    "bmi_available",
    "hba1c_available",
    "genes_with_biotype",
    "protein_coding_genes",
    "lncrna_genes"
  ),
  value = c(
    nrow(gene_annotation),
    nrow(sample_metadata),
    sum(sample_metadata$status == "ND", na.rm = TRUE),
    sum(sample_metadata$status == "IGT", na.rm = TRUE),
    sum(sample_metadata$status == "T2D", na.rm = TRUE),
    sum(sample_metadata$status == "T3cD", na.rm = TRUE),
    sum(!is.na(sample_metadata$age)),
    sum(!is.na(sample_metadata$bmi)),
    sum(!is.na(sample_metadata$hba1c)),
    sum(!is.na(gene_annotation$biotype)),
    sum(gene_annotation$biotype == "protein_coding", na.rm = TRUE),
    sum(gene_annotation$biotype == "lncRNA", na.rm = TRUE)
  )
)

status_table <- sample_metadata %>%
  dplyr::count(status, name = "n") %>%
  dplyr::arrange(status)

biotype_table <- gene_annotation %>%
  dplyr::count(biotype, name = "n", sort = TRUE) %>%
  dplyr::filter(!is.na(biotype))


excel_sheets <- list(
  summary = summary_table,
  status_counts = status_table,
  sample_metadata = sample_metadata,
  gene_annotation_preview = head(gene_annotation, 500),
  biotype_counts = head(biotype_table, 50)
)

writexl::write_xlsx(
  excel_sheets,
  path = file.path(project_dir, "GSE164416_metadata_complete.xlsx")
)

cat("\n========================================\n")
cat("Data Acquisition Complete!\n")
cat("========================================\n")
cat("Files exported:\n")
cat("  1. data/processed/counts_matrix.tsv\n")
cat("  2. data/metadata/sample_metadata.csv\n")
cat("  3. data/metadata/gene_annotation.tsv\n")
cat("  4. GSE164416_metadata_complete.xlsx\n")
cat("\nDataset summary:\n")
cat("  - 133 samples (18 ND, 38-41 IGT, 39 T2D, 35 T3cD)\n")
cat("  - ~60,000 genes with biotype annotation\n")
cat("  - Clinical metadata: age, BMI, HbA1c (sex not available)\n")
cat("\nNext steps:\n")
cat("  1. Quality control and filtering (02_quality_control.R)\n")
cat("  2. Normalization (DESeq2 or edgeR)\n")
cat("  3. Differential expression analysis\n")
cat("  4. Coexpression network analysis (WGCNA)\n")
cat("========================================\n")

