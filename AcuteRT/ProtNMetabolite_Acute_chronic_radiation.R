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
    tp1, tp2,                        # e.g. "B0", "F3"
    file_prefix,                     # "protein" or "metabolite"
    sample_col = "Sample_ID",
    time_col   = "Time",
    arm_col    = "Arm",              # NEW
    arm_include = "NIF",             # NEW (can be NULL to include all)
    outdir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/LimmaResults/ExtraSig/"
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
  req_cols <- c(sample_col, time_col)
  if (!all(req_cols %in% colnames(metaData))) {
    stop("metaData must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  metaSubset <- metaData %>%
    dplyr::filter(.data[[time_col]] %in% c(tp1, tp2))
  
  # Keep only NIF (configurable). If arm_include is NULL, skip this.
  if (!is.null(arm_include)) {
    if (!arm_col %in% colnames(metaSubset)) {
      stop("Requested arm filtering but '", arm_col, "' not found in metaData.")
    }
    metaSubset <- metaSubset %>% dplyr::filter(.data[[arm_col]] %in% arm_include)
  }
  
  # Patient_ID (if you use within-patient blocking)
  metaSubset <- metaSubset %>%
    dplyr::mutate(Patient_ID = sub("_.*", "", .data[[sample_col]]))
  
  # Align to expression data
  common <- intersect(colnames(data_num), metaSubset[[sample_col]])
  if (length(common) < 2) stop("Not enough overlapping samples after filtering (time + arm).")
  metaSubset <- metaSubset[match(common, metaSubset[[sample_col]]), , drop = FALSE]
  data_subset <- data_num[, common, drop = FALSE]
  
  # Factor with baseline as reference
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
  
  # Optional significance filter (your current choice)
  res <- res[res$P.Value <= 0.05 & abs(res$logFC) >= 1, , drop = FALSE]
  
  # --- Merge with molecule metadata ---
  res$RowID <- rownames(res)
  merged_res <- if (tolower(file_prefix) == "protein") {
    dplyr::left_join(tibble::as_tibble(res, rownames = NULL), molMetadata, by = c("RowID" = "SeqID"))
  } else if (tolower(file_prefix) == "metabolite") {
    dplyr::left_join(tibble::as_tibble(res, rownames = NULL), molMetadata, by = c("RowID" = "Index"))
  } else {
    stop("file_prefix must be 'protein' or 'metabolite'.")
  }
  
  # --- Save CSV (add arm tag if filtered) ---
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  arm_tag <- if (is.null(arm_include)) "ALL" else paste(arm_include, collapse = "-")
  outfile <- file.path(outdir, sprintf("%s_%s_vs_%s_%s.csv", file_prefix, tp2, tp1, arm_tag))
  utils::write.csv(merged_res, outfile, row.names = FALSE)
  message("Saved: ", outfile)
  
  merged_res
}




plot_feature_trends <- function(
    feature_names,
    output_dir,
    expr_data,
    metaData,
    sample_col  = "Sample_ID",
    time_col    = "Time",
    arm_col     = "Arm",          # NEW
    arm_include = "NIF",          # NEW (set NULL to include all arms)
    valid_times = c("B0","D5","F3"),
    transform_y = c("none","log2p1"),
    molMetadata = NULL,
    id_col      = NULL,
    title_col   = NULL
){
  transform_y <- match.arg(transform_y)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Folder created: ", output_dir)
  } else {
    message("Folder exists: ", output_dir)
  }
  
  # keep only numeric sample columns
  num_cols <- vapply(expr_data, is.numeric, logical(1))
  X <- as.matrix(expr_data[, num_cols, drop = FALSE])
  
  # optional transform
  if (transform_y == "log2p1") {
    X <- log2(X + 1)
    ylab <- "Log2(+1) abundance"
  } else {
    ylab <- "Abundance"
  }
  
  # basic checks
  stopifnot(all(c(sample_col, time_col) %in% colnames(metaData)))
  
  # ---- Arm filter (before alignment) ----
  meta_filt <- metaData
  if (!is.null(arm_include)) {
    if (!arm_col %in% colnames(meta_filt)) {
      stop("Requested arm filtering but '", arm_col, "' not found in metaData.")
    }
    meta_filt <- dplyr::filter(meta_filt, .data[[arm_col]] %in% arm_include)
  }
  
  # align to expression data
  common <- intersect(colnames(X), meta_filt[[sample_col]])
  if (length(common) < 2) {
    stop("Not enough overlapping samples between expr_data and metaData after arm filtering.")
  }
  meta_sub <- meta_filt[match(common, meta_filt[[sample_col]]), , drop = FALSE]
  X <- X[, common, drop = FALSE]
  
  # filter by valid times
  keep_time <- meta_sub[[time_col]] %in% valid_times
  if (!any(keep_time)) stop("No samples in requested valid_times: ", paste(valid_times, collapse=", "))
  meta_sub <- meta_sub[keep_time, , drop = FALSE]
  X <- X[, keep_time, drop = FALSE]
  meta_sub[[time_col]] <- factor(meta_sub[[time_col]], levels = valid_times)
  
  # ---- build a lookup from molMetadata if provided ----
  title_lookup <- NULL
  if (!is.null(molMetadata) && !is.null(id_col) && !is.null(title_col)) {
    if (!all(c(id_col, title_col) %in% colnames(molMetadata))) {
      stop("molMetadata must contain columns: ", paste(c(id_col, title_col), collapse = ", "))
    }
    key   <- as.character(molMetadata[[id_col]])
    value <- as.character(molMetadata[[title_col]])
    title_lookup <- stats::setNames(value[!duplicated(key)], key[!duplicated(key)])
  }
  
  # iterate features
  feature_names <- as.character(feature_names)
  for (feat in feature_names) {
    if (!feat %in% rownames(X)) {
      warning("Feature not found in expr_data: ", feat); next
    }
    
    pretty_title <- if (!is.null(title_lookup) && feat %in% names(title_lookup)) title_lookup[[feat]] else feat
    
    vals <- as.numeric(X[feat, , drop = TRUE])
    df <- data.frame(Value = vals, stringsAsFactors = FALSE, check.names = FALSE)
    df[[sample_col]] <- colnames(X)
    df <- df[, c(sample_col, "Value")]
    
    df <- dplyr::left_join(df, meta_sub, by = setNames(sample_col, sample_col))
    df <- dplyr::filter(df, .data[[time_col]] %in% valid_times)
    if (nrow(df) == 0) { warning("No data for feature ", feat, " at requested timepoints."); next }
    df[[time_col]] <- factor(df[[time_col]], levels = valid_times)
    
    df_sum <- df %>%
      dplyr::group_by(.data[[time_col]]) %>%
      dplyr::summarise(
        meanValue = mean(Value, na.rm = TRUE),
        sdValue   = sd(Value, na.rm = TRUE),
        n         = dplyr::n(),
        .groups   = "drop"
      )
    
    p <- ggplot(df, aes(x = .data[[time_col]], y = Value)) +
      geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
      geom_point(data = df_sum, aes(x = .data[[time_col]], y = meanValue), size = 3, shape = 18, inherit.aes = FALSE) +
      geom_line (data = df_sum, aes(x = .data[[time_col]], y = meanValue, group = 1), linewidth = 0.6, inherit.aes = FALSE) +
      geom_errorbar(data = df_sum, aes(x = .data[[time_col]], ymin = meanValue - sdValue, ymax = meanValue + sdValue),
                    width = 0.15, inherit.aes = FALSE) +
      labs(title = pretty_title, x = "Time", y = ylab, subtitle = feat) +
      theme_bw()
    
    file_suffix <- if (is.null(arm_include)) "ALL" else paste(arm_include, collapse = "-")
    file_out <- file.path(output_dir, paste0(feat, "_trend_", file_suffix, ".png"))
    ggsave(file_out, plot = p, width = 8, height = 6, dpi = 300)
    message("Saved plot: ", file_out)
  }
}

# =============================================================================
# Stats analysis

# Focused Analysis
prot_NIF <- run_limma_timepoints(
  data        = protData,
  metaData    = prot_metaData,
  molMetadata = protein_mol_metadata,
  tp1         = "B0",
  tp2         = "D1",
  file_prefix = "protein",
  arm_include = "NIF",
  outdir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/LimmaResults/ExtraSig_NIF/"
)

metl_results <- run_limma_timepoints(
  data        = metlData,
  metaData    = metl_metaData,
  molMetadata = metabolite_mol_meta_data,
  tp1         = "B0",
  tp2         = "D1",
  file_prefix = "metabolite",
  arm_include = "NIF",
  outdir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/LimmaResults/ExtraSig_NIF/"
)


# =============================================================================
# Figures

# Glycerophospholipids
plot_feature_trends(
  feature_names = c("MW0055322"),
  output_dir    = out_figs_metl,
  expr_data     = metlData,
  metaData      = metl_metaData,
  arm_include   = "IF",                 # <- filter to NIF
  valid_times   = c("B0","D1","D5","F3", "F6"),
  transform_y   = "log2p1",
  molMetadata   = metabolite_mol_meta_data,
  id_col        = "Index",
  title_col     = "Compounds"
)

plot_feature_trends(
  feature_names = c("MW0060249","MW0059344", "MEDN1278", "MEDP1322"),
  output_dir    = out_figs_metl,
  expr_data     = metlData,
  metaData      = metl_metaData,
  arm_include   = "IF",                 # <- filter to NIF
  valid_times   = c("B0","D1","D5", "F3", "F6"),
  transform_y   = "log2p1",
  molMetadata   = metabolite_mol_meta_data,
  id_col        = "Index",
  title_col     = "Compounds"
)


# Lysophosphatidylcholines
plot_feature_trends(
  feature_names = c("MEDN1278","MEDP1322"),
  output_dir    = out_figs_metl,
  expr_data     = metlData,
  metaData      = metl_metaData,
  arm_include   = "NIF",                 # <- filter to NIF
  valid_times   = c("B0","D1","F3"),
  transform_y   = "log2p1",
  molMetadata   = metabolite_mol_meta_data,
  id_col        = "Index",
  title_col     = "Compounds"
)

#Prostaglandin 

plot_feature_trends(
  feature_names = c("MEDN1430","MW0015050"),
  output_dir    = out_figs_metl,
  expr_data     = metlData,
  metaData      = metl_metaData,
  arm_include   = "NIF",                 # <- filter to NIF
  valid_times   = c("B0","D1","D5", "F3", "F6"),
  transform_y   = "log2p1",
  molMetadata   = metabolite_mol_meta_data,
  id_col        = "Index",
  title_col     = "Compounds"
)

# =============================================================================
# Legacy

#########################################

# Protein example
prot_results <- run_limma_timepoints(
  data        = protData,
  metaData    = prot_metaData,
  molMetadata = protein_mol_metadata,
  tp1         = "B0",
  tp2         = "D1",
  file_prefix = "protein"
)

prot_results <- run_limma_timepoints(
  data        = protData,
  metaData    = prot_metaData,
  molMetadata = protein_mol_metadata,
  tp1         = "D1",
  tp2         = "F3",
  file_prefix = "protein"
)

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
  tp2         = "D1",
  file_prefix = "metabolite"
)

metl_results <- run_limma_timepoints(
  data        = metlData,
  metaData    = metl_metaData,
  molMetadata = metabolite_mol_meta_data,
  tp1         = "D1",
  tp2         = "F3",
  file_prefix = "metabolite"
)

metl_results <- run_limma_timepoints(
  data        = metlData,
  metaData    = metl_metaData,
  molMetadata = metabolite_mol_meta_data,
  tp1         = "B0",
  tp2         = "F3",
  file_prefix = "metabolite"
)


#---------------------------- Ploting ---------------------------------------

out_figs_prot <- "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/figs/prot"
out_figs_metl <- "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/figs/metl"

plot_feature_trends(
  feature_names = "MW0105369",
  output_dir    = out_figs_metl,
  expr_data     = metlData,
  metaData      = metl_metaData,
  valid_times   = c("B0","D1","F3"),
  transform_y   = "log2p1"
)

plot_feature_trends(
  feature_names = "seq.6895.1",
  output_dir    = out_figs_prot,
  expr_data     = protData,
  metaData      = prot_metaData,
  valid_times   = c("B0","D5","F3", "F6"),
  transform_y   = "log2p1"
)

plot_feature_trends(
  feature_names = "seq.21796.43",
  output_dir    = out_figs_prot,
  expr_data     = protData,
  metaData      = prot_metaData,
  valid_times   = c("B0","D1","F3"),
  transform_y   = "log2p1",
  molMetadata   = protein_mol_metadata,
  id_col        = "SeqID",
  title_col     = "Description"
)

plot_feature_trends(
  feature_names = "seq.21796.43",
  output_dir    = out_figs_prot,
  expr_data     = protData,
  metaData      = prot_metaData,
  valid_times   = c("B0","D1","F3"),
  transform_y   = "none",
    molMetadata   = protein_mol_metadata,
  id_col        = "SeqID",
  title_col     = "Description"
)
