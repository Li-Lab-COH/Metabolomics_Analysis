# =============================================================================


# Load required libraries
library(limma)
library(dplyr)
library(readxl)
library(stringr)
library(caret)
library(matrixStats)
library(ggplot2)

# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12

# =============================================================================
# Loading data

# Loading Protein data
protData <- read.csv("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/ResultsTables_Filtered_1.csv", row.names = "SeqID")
prot_metaData <- read_excel("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/meta_data_proteomics 2.xlsx")
protein_mol_metadata <- read.csv("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/protein_metadata.csv")

#Loading Metabolite data
# metlData <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_data_74.csv", row.names = "ID")
metlData <- read.csv(
  "~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_data_74.csv",
  row.names = "ID",
  check.names = FALSE
)
metl_metaData <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/meta_data_74.csv")
metabolite_mol_meta_data <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_metadata.csv")
rownames(metlData) <- sub("\\|.*", "", rownames(metlData))


# =============================================================================
# Function for testing
run_limma_timepoints <- function(
    data,
    metaData,
    molMetadata,
    tp1, tp2,                  # e.g. "B0", "F3"
    file_prefix,               # "protein" or "metabolite"
    sample_col = "Sample_ID",
    time_col   = "Time",
    outdir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT"
){
  # --- Process data ---
  num_cols <- vapply(data, is.numeric, logical(1))
  if (!any(num_cols)) stop("No numeric columns in 'data'.")
  data_num <- as.matrix(data[, num_cols, drop = FALSE])
  data_num <- log2(data_num + 1)
  
  vars <- matrixStats::rowVars(data_num, na.rm = TRUE)
  var_cut <- quantile(vars, 0.15, na.rm = TRUE)
  data_num <- data_num[vars > var_cut, , drop = FALSE]
  
  # --- Subset metadata ---
  if (!all(c(sample_col, time_col) %in% colnames(metaData))) {
    stop("metaData must contain columns: ", paste(c(sample_col, time_col), collapse = ", "))
  }
  metaSubset <- metaData %>%
    dplyr::filter(.data[[time_col]] %in% c(tp1, tp2)) %>%
    dplyr::mutate(Patient_ID = sub("_.*", "", .data[[sample_col]]))
  
  common <- intersect(colnames(data_num), metaSubset[[sample_col]])
  if (length(common) < 2) stop("Not enough overlapping samples.")
  metaSubset <- metaSubset[match(common, metaSubset[[sample_col]]), , drop = FALSE]
  data_subset <- data_num[, common, drop = FALSE]
  
  metaSubset[[time_col]] <- factor(metaSubset[[time_col]], levels = c(tp1, tp2))
  
  # --- Limma model ---
  design <- model.matrix(stats::as.formula(paste0("~ ", time_col)), data = metaSubset)
  corfit <- limma::duplicateCorrelation(data_subset, design, block = metaSubset$Patient_ID)
  fit <- limma::lmFit(data_subset, design, block = metaSubset$Patient_ID, correlation = corfit$consensus)
  fit <- limma::eBayes(fit)
  
  coef_name <- paste0(time_col, tp2)
  if (!coef_name %in% colnames(fit$coefficients)) {
    stop("Could not find coefficient '", coef_name, "'.")
  }
  
  res <- limma::topTable(fit, coef = coef_name, number = Inf)
  res <- res[res$P.Value <= 0.05, , drop = FALSE]
  
  # --- Merge with molecule metadata ---
  res$RowID <- rownames(res)
  if (tolower(file_prefix) == "protein") {
    # Match RowID to SeqID in protein metadata
    merged_res <- dplyr::left_join(
      tibble::as_tibble(res, rownames = NULL),
      molMetadata,
      by = c("RowID" = "SeqID")
    )
  } else if (tolower(file_prefix) == "metabolite") {
    # Match RowID to Index in metabolite metadata
    merged_res <- dplyr::left_join(
      tibble::as_tibble(res, rownames = NULL),
      molMetadata,
      by = c("RowID" = "Index")
    )
  } else {
    stop("file_prefix must be 'protein' or 'metabolite'.")
  }
  
  # --- Save CSV ---
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outdir, sprintf("%s_%s_vs_%s.csv", file_prefix, tp2, tp1))
  utils::write.csv(merged_res, outfile, row.names = FALSE)
  message("Saved: ", outfile)
  
  return(merged_res)
}

# Protein example
prot_results <- run_limma_timepoints(
  data        = protData,
  metaData    = prot_metaData,
  molMetadata = protein_mol_metadata,
  tp1         = "B0",
  tp2         = "F3",
  file_prefix = "protein"
)

# Metabolite example
metl_results <- run_limma_timepoints(
  data        = metlData,
  metaData    = metl_metaData,
  molMetadata = metabolite_mol_meta_data,
  tp1         = "B0",
  tp2         = "F3",
  file_prefix = "metabolite"
)

