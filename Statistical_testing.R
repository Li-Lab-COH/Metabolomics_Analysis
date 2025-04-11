# ================================
# 0. Setup & Directory Creation
# ================================
"~/Roselab/Metabolite/results/"
# Function to create directories if they do not exist
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}

# Define all output directories
dirs <- c(
  "~/Roselab/Metabolite/results/difference/limma/",
  "~/Roselab/Metabolite/results/difference/limma/withinArm/",
  "~/Roselab/Metabolite/results/difference/limma/intraArm/",
  "~/Roselab/Metabolite/results/difference/t_test/",
  "~/Roselab/Metabolite/results/difference/t_test/withinArm/",
  "~/Roselab/Metabolite/results/difference/t_test/intraArm/"
)
sapply(dirs, ensure_dir)

# ================================
# 1. Data Loading & Pre-processing
# ================================
library(limma)
library(dplyr)
# (If not using metaboanalystR for the t.test analysis, you may remove it)
# library(MetaboAnalystR)

# Load data files
metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID")
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")
metabolite_meta_data <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_metadata.csv")

# Log2 transform (adding 1 to avoid log(0))
# metlData_log2 <- log2(metlData + 1)
metlData_log2_filtered <- log2(metlData + 1)

# # Filter out metabolites with low variance
# row_vars <- apply(metlData_log2, 1, var)
# variance_threshold <- 0.5
# metlData_log2_filtered <- metlData_log2[row_vars >= variance_threshold, ]

cat("Before filtering:", nrow(metlData_log2), "metabolites\n")
cat("After filtering:", nrow(metlData_log2_filtered), "metabolites\n")

# Define time points and subset metadata
timepoints_of_interest <- c("B0", "D5", "F3", "F6")
metaSubset <- metaData %>%
  filter(Time %in% timepoints_of_interest, Arm %in% c("IF", "NIF")) %>% 
  mutate(Patient_ID = gsub("_.*", "", Sample_ID))
metaSubset$Time <- factor(metaSubset$Time, levels = timepoints_of_interest)
metaSubset$Arm <- factor(metaSubset$Arm, levels = c("NIF", "IF"))

# ------------------------------
# Function to process results:
#   - Split the rownames by "|" into Index and Compound
#   - Merge with metabolite metadata (by column "Index")
# ------------------------------
process_results <- function(results_df, metabolite_meta_data) {
  # Split row names (assumed to be "Index|Compound")
  split_names <- strsplit(rownames(results_df), "\\|")
  index <- sapply(split_names, `[`, 1)
  compound <- sapply(split_names, function(x) paste(x[-1], collapse="|"))
  results_df <- cbind(Index = index, Compound = compound, results_df)
  # Merge with metabolite metadata (assuming metabolite_meta_data has a column named "Index")
  merged_df <- merge(metabolite_meta_data, results_df, by = "Index", all.y = TRUE)
  return(merged_df)
}

# ================================
# 2. LIMMA Analysis
# ================================

## (A) IntraArm Comparisons (Between groups at the same time point)
for(tp in timepoints_of_interest){
  message("Running limma intraArm analysis for time point: ", tp)
  
  # Subset metadata and corresponding metabolite data
  meta_tp <- metaSubset[metaSubset$Time == tp, ]
  metlData_tp <- metlData_log2_filtered[, meta_tp$Sample_ID, drop = FALSE]
  
  # Design matrix with Arm as the only variable (NIF as reference)
  design_tp <- model.matrix(~ Arm, data = meta_tp)
  
  # Fit the limma model
  fit_tp <- lmFit(metlData_tp, design_tp)
  fit_tp <- eBayes(fit_tp)
  
  # Contrast: difference IF vs NIF (coefficient "ArmIF" if NIF is the reference)
  results_tp <- topTable(fit_tp, coef = "ArmIF", number = Inf)
  
  # Process results: split row names and merge with metabolite metadata
  processed_tp <- process_results(results_tp, metabolite_meta_data)
  
  # Save final merged output
  outfile <- file.path("~/Roselab/Metabolite/results/difference/limma/intraArm",
                       paste0("limma_intraArm_", tp, "_IF_vs_NIF.csv"))
  write.csv(processed_tp, file = outfile, row.names = FALSE)
  message("Results saved: ", outfile)
}

## (B) WithinArm Comparisons (Within each arm: later time points vs baseline B0)
arms <- c("NIF", "IF")
for(arm in arms){
  message("Running limma withinArm analysis for arm: ", arm)
  meta_arm <- metaSubset[metaSubset$Arm == arm, ]
  
  for(tp in setdiff(timepoints_of_interest, "B0")){
    message("  Comparing B0 vs ", tp, " for arm ", arm)
    
    # Subset for baseline and the given time point
    meta_arm_tp <- meta_arm[meta_arm$Time %in% c("B0", tp), ]
    metlData_arm_tp <- metlData_log2_filtered[, meta_arm_tp$Sample_ID, drop = FALSE]
    
    # Order Time factor with baseline as reference
    meta_arm_tp$Time <- factor(meta_arm_tp$Time, levels = c("B0", tp))
    
    # Build design matrix
    design_arm <- model.matrix(~ Time, data = meta_arm_tp)
    
    # Estimate within-patient correlation
    corfit_arm <- duplicateCorrelation(metlData_arm_tp, design_arm, block = meta_arm_tp$Patient_ID)
    
    # Fit the model using the estimated correlation
    fit_arm <- lmFit(metlData_arm_tp, design_arm, block = meta_arm_tp$Patient_ID, 
                     correlation = corfit_arm$consensus)
    fit_arm <- eBayes(fit_arm)
    
    # Coefficient name, e.g., "TimeD5", "TimeF3", etc.
    coef_name <- paste0("Time", tp)
    results_arm <- topTable(fit_arm, coef = coef_name, number = Inf)
    
    # Process results with metadata
    processed_arm <- process_results(results_arm, metabolite_meta_data)
    
    # Save output
    outfile <- file.path("~/Roselab/Metabolite/results/difference/limma/withinArm",
                         paste0("limma_withinArm_", arm, "_B0_vs_", tp, ".csv"))
    write.csv(processed_arm, file = outfile, row.names = FALSE)
    message("  Results saved: ", outfile)
  }
}

# ================================
# 3. T.TEST Analysis
# ================================
# For the t.test analyses, we perform similar comparisons as above.
# A helper function is defined to run t.tests for each metabolite (row) in a data matrix.
#
# For intraArm comparisons, an unpaired t.test is performed.
# For withinArm comparisons, we perform paired t.test using only subjects having paired samples.

# Helper function to perform t.test for one metabolite given a numeric vector x and grouping vector group
run_ttest <- function(x, group, paired = FALSE, timePoint) {
  # Ensure x is numeric; group is a factor with two levels.
  # Use t.test with specified pairing if needed.
  
  # Split values by group
  group1 <- x[group == levels(group)[1]]
  group2 <- x[group == levels(group)[2]]

  # Debug print
  cat("\n----- Running t.test - time", tp,  "-----\n")
  cat("Group 1 (", levels(group)[1], ") values:\n")
  print(group1)
  cat("Group 2 (", levels(group)[2], ") values:\n")
  print(group2)

  # Print corresponding sample names
  cat("Sample IDs for Group 1:\n")
  print(names(group1))
  cat("Sample IDs for Group 2:\n")
  print(names(group2))


  if(paired){
    # For paired t-test, assume that the order in group must match for pairs.
    tt <- t.test(x[group == levels(group)[1]], x[group == levels(group)[2]], paired = TRUE)
  } else {
    tt <- t.test(x ~ group)
  }
  # Return a named vector with relevant statistics
  c(t_stat = tt$statistic, p_value = tt$p.value, 
    mean_group1 = ifelse("mean in group 1" %in% names(tt$estimate), tt$estimate[1], NA),
    mean_group2 = ifelse("mean in group 2" %in% names(tt$estimate), tt$estimate[2], NA))
}

## (A) T.test IntraArm Comparisons (IF vs NIF at each time point)
intra_ttest_results <- list()
timepoints_of_interest_ttest <- c("D5", "F3", "F6")
for(tp in timepoints_of_interest_ttest){
  message("Running t.test intraArm analysis for time point: ", tp)
  
  meta_tp <- metaSubset[metaSubset$Time == tp, ]
  print(meta_tp)
  metlData_tp <- metlData_log2_filtered[, meta_tp$Sample_ID, drop = FALSE]
  
  # Group vector from metadata
  group_vec <- factor(as.character(meta_tp$Arm), levels = c("NIF", "IF"))
  print("Here is the group_vec")
  print(group_vec)
  print(head(metlData_tp))
  
  print("Removing metabolites with 0 variance")
  # Filter out metabolites (rows) with zero variance in either group
  non_constant_metlData_tp <- metlData_tp[apply(metlData_tp, 1, function(x) {
    group1 <- x[group_vec == levels(group_vec)[1]]
    group2 <- x[group_vec == levels(group_vec)[2]]
    var(group1, na.rm = TRUE) > 0 && var(group2, na.rm = TRUE) > 0
  }), ]
  
  # Apply t.test for each metabolite (each row)
  ttest_out <- t(apply(metlData_tp, 1, run_ttest, group = group_vec, paired = FALSE, timePoint = tp))
  ttest_out <- as.data.frame(ttest_out)
  
  # For consistency, add row names as a column and process them:
  ttest_out$RowID <- rownames(ttest_out)
  rownames(ttest_out) <- NULL
  
  # Split the RowID into Index and Compound
  split_names <- strsplit(ttest_out$RowID, "\\|")
  ttest_out$Index <- sapply(split_names, `[`, 1)
  ttest_out$Compound <- sapply(split_names, function(x) paste(x[-1], collapse="|"))
  ttest_out$RowID <- NULL
  
  # Merge with metabolite metadata
  merged_ttest <- merge(metabolite_meta_data, ttest_out, by = "Index", all.y = TRUE)
  
  # Save results file
  outfile <- file.path("~/Roselab/Metabolite/results/difference/t_test/intraArm",
                       paste0("t_test_intraArm_", tp, "_IF_vs_NIF.csv"))
  write.csv(merged_ttest, file = outfile, row.names = FALSE)
  message("T.test results saved: ", outfile)
}

## (B) T.test WithinArm Comparisons (For each arm, compare B0 vs later time point via paired t.test)
arms <- c("NIF", "IF")
for(arm in arms){
  message("Running t.test withinArm analysis for arm: ", arm)
  meta_arm <- metaSubset[metaSubset$Arm == arm, ]
  
  for(tp in setdiff(timepoints_of_interest, "B0")){
    message("  Comparing B0 vs ", tp, " for arm ", arm)
    
    # Subset to samples with either baseline (B0) or current tp
    meta_arm_tp <- meta_arm[meta_arm$Time %in% c("B0", tp), ]
    # Identify patients that have samples at both timepoints
    common_patients <- intersect(
      meta_arm_tp$Patient_ID[meta_arm_tp$Time == "B0"],
      meta_arm_tp$Patient_ID[meta_arm_tp$Time == tp]
    )
    meta_arm_tp <- meta_arm_tp[meta_arm_tp$Patient_ID %in% common_patients, ]
    
    # Order the factor so that B0 is level1 and tp is level2
    meta_arm_tp$Time <- factor(meta_arm_tp$Time, levels = c("B0", tp))
    
    metData_arm_tp <- metlData_log2_filtered[, meta_arm_tp$Sample_ID, drop = FALSE]
    
    # For a paired t.test, we need to match the samples by Patient_ID.
    # Create separate data matrices for B0 and the current timepoint.
    B0_samples <- meta_arm_tp$Sample_ID[meta_arm_tp$Time == "B0"]
    tp_samples <- meta_arm_tp$Sample_ID[meta_arm_tp$Time == tp]
    
    # Ensure proper ordering by Patient_ID
    B0_meta <- meta_arm_tp[meta_arm_tp$Time == "B0", ]
    tp_meta <- meta_arm_tp[meta_arm_tp$Time == tp, ]
    B0_meta <- B0_meta[order(B0_meta$Patient_ID), ]
    tp_meta <- tp_meta[order(tp_meta$Patient_ID), ]
    B0_samples <- B0_meta$Sample_ID
    tp_samples <- tp_meta$Sample_ID
    
    # Initialize matrix to hold t.test results for each metabolite
    ttest_res <- t(apply(metData_arm_tp, 1, function(x) {
      # Get paired values (ensure that both samples exist for a given metabolite for each patient)
      x_B0 <- as.numeric(x[B0_samples])
      x_tp <- as.numeric(x[tp_samples])
      # perform paired t-test
      tt <- t.test(x_tp, x_B0, paired = TRUE)
      c(t_stat = tt$statistic, p_value = tt$p.value, 
        mean_B0 = mean(x_B0, na.rm = TRUE), mean_tp = mean(x_tp, na.rm = TRUE))
    }))
    
    ttest_res <- as.data.frame(ttest_res)
    ttest_res$RowID <- rownames(ttest_res)
    rownames(ttest_res) <- NULL
    
    # Split RowID into Index and Compound
    split_names <- strsplit(ttest_res$RowID, "\\|")
    ttest_res$Index <- sapply(split_names, `[`, 1)
    ttest_res$Compound <- sapply(split_names, function(x) paste(x[-1], collapse = "|"))
    ttest_res$RowID <- NULL
    
    # Merge with metabolite metadata
    merged_ttest <- merge(metabolite_meta_data, ttest_res, by = "Index", all.y = TRUE)
    
    # Save output
    outfile <- file.path("~/Roselab/Metabolite/results/difference/t_test/withinArm",
                         paste0("t_test_withinArm_", arm, "_B0_vs_", tp, ".csv"))
    write.csv(merged_ttest, file = outfile, row.names = FALSE)
    message("  T.test results saved: ", outfile)
  }
}

# ================================
# End of Script
# ================================












# ---------------------------
# Helper function to perform t.test for one metabolite given a numeric vector x and grouping vector group
# ---------------------------
run_ttest <- function(x, group, paired = FALSE) {
  # Ensure x is numeric; group is a factor with two levels.
  # Use t.test with specified pairing if needed.
  if(paired){
    # For paired t-test, assume that the order in group must match for pairs.
    tt <- t.test(x[group == levels(group)[1]], x[group == levels(group)[2]], paired = TRUE)
  } else {
    tt <- t.test(x ~ group)
  }
  # Return a named vector with relevant statistics
  c(t_stat = tt$statistic, p_value = tt$p.value, 
    mean_group1 = ifelse("mean in group 1" %in% names(tt$estimate), tt$estimate[1], NA),
    mean_group2 = ifelse("mean in group 2" %in% names(tt$estimate), tt$estimate[2], NA))
}

# ---------------------------
# (A) T.test IntraArm Comparisons (IF vs NIF at each time point)
# ---------------------------
intra_ttest_results <- list()
for(tp in timepoints_of_interest){
  message("Running t.test intraArm analysis for time point: ", tp)
  
  meta_tp <- metaSubset[metaSubset$Time == tp, ]
  metlData_tp <- metlData_log2_filtered[, meta_tp$Sample_ID, drop = FALSE]
  
  # Group vector from metadata (NIF vs IF)
  group_vec <- factor(as.character(meta_tp$Arm), levels = c("NIF", "IF"))
  
  # Apply t.test for each metabolite (each row)
  ttest_out <- t(apply(metlData_tp, 1, run_ttest, group = group_vec, paired = FALSE))
  ttest_out <- as.data.frame(ttest_out)
  
  # Explicitly assign meaningful column names for intraArm results:
  colnames(ttest_out)[1:4] <- c("t_stat", "p_value", "mean_NIF", "mean_IF")
  
  # For consistency, add row names as a column and process them:
  ttest_out$RowID <- rownames(ttest_out)
  rownames(ttest_out) <- NULL
  
  # Split the RowID into Index and Compound (the original row names contain "Index|Compound")
  split_names <- strsplit(ttest_out$RowID, "\\|")
  ttest_out$Index <- sapply(split_names, `[`, 1)
  ttest_out$Compound <- sapply(split_names, function(x) paste(x[-1], collapse = "|"))
  ttest_out$RowID <- NULL
  
  # Merge with metabolite metadata
  merged_ttest <- merge(metabolite_meta_data, ttest_out, by = "Index", all.y = TRUE)
  
  # Save results file
  outfile <- file.path("~/Roselab/Metabolite/results/difference/t_test/intraArm",
                       paste0("t_test_intraArm_", tp, "_IF_vs_NIF.csv"))
  write.csv(merged_ttest, file = outfile, row.names = FALSE)
  message("T.test intraArm results saved: ", outfile)
}

# ---------------------------
# (B) T.test WithinArm Comparisons (For each arm, compare B0 vs later time point via paired t.test)
# ---------------------------
arms <- c("NIF", "IF")
for(arm in arms){
  message("Running t.test withinArm analysis for arm: ", arm)
  meta_arm <- metaSubset[metaSubset$Arm == arm, ]
  
  for(tp in setdiff(timepoints_of_interest, "B0")){
    message("  Comparing B0 vs ", tp, " for arm ", arm)
    
    # Subset to samples with either baseline (B0) or current timepoint
    meta_arm_tp <- meta_arm[meta_arm$Time %in% c("B0", tp), ]
    # Identify patients that have samples at both timepoints
    common_patients <- intersect(
      meta_arm_tp$Patient_ID[meta_arm_tp$Time == "B0"],
      meta_arm_tp$Patient_ID[meta_arm_tp$Time == tp]
    )
    meta_arm_tp <- meta_arm_tp[meta_arm_tp$Patient_ID %in% common_patients, ]
    
    # Order the factor so that B0 is level1 and tp is level2
    meta_arm_tp$Time <- factor(meta_arm_tp$Time, levels = c("B0", tp))
    
    metData_arm_tp <- metlData_log2_filtered[, meta_arm_tp$Sample_ID, drop = FALSE]
    
    # For a paired t.test, we need to match the samples by Patient_ID.
    # Create separate data matrices for B0 and the current timepoint.
    B0_samples <- meta_arm_tp$Sample_ID[meta_arm_tp$Time == "B0"]
    tp_samples <- meta_arm_tp$Sample_ID[meta_arm_tp$Time == tp]
    
    # Ensure proper ordering by Patient_ID
    B0_meta <- meta_arm_tp[meta_arm_tp$Time == "B0", ]
    tp_meta <- meta_arm_tp[meta_arm_tp$Time == tp, ]
    B0_meta <- B0_meta[order(B0_meta$Patient_ID), ]
    tp_meta <- tp_meta[order(tp_meta$Patient_ID), ]
    B0_samples <- B0_meta$Sample_ID
    tp_samples <- tp_meta$Sample_ID
    
    # Initialize matrix to hold t.test results for each metabolite, performed as paired t-test.
    ttest_res <- t(apply(metData_arm_tp, 1, function(x) {
      x_B0 <- as.numeric(x[B0_samples])
      x_tp <- as.numeric(x[tp_samples])
      tt <- t.test(x_tp, x_B0, paired = TRUE)
      c(t_stat = tt$statistic, p_value = tt$p.value, 
        mean_B0 = mean(x_B0, na.rm = TRUE), mean_tp = mean(x_tp, na.rm = TRUE))
    }))
    
    ttest_res <- as.data.frame(ttest_res)
    # Explicitly assign column names for withinArm results:
    colnames(ttest_res)[1:4] <- c("t_stat", "p_value", "mean_B0", "mean_tp")
    
    ttest_res$RowID <- rownames(ttest_res)
    rownames(ttest_res) <- NULL
    
    # Split RowID into Index and Compound
    split_names <- strsplit(ttest_res$RowID, "\\|")
    ttest_res$Index <- sapply(split_names, `[`, 1)
    ttest_res$Compound <- sapply(split_names, function(x) paste(x[-1], collapse = "|"))
    ttest_res$RowID <- NULL
    
    # Merge with metabolite metadata
    merged_ttest <- merge(metabolite_meta_data, ttest_res, by = "Index", all.y = TRUE)
    
    # Save output
    outfile <- file.path("~/Roselab/Metabolite/results/difference/t_test/withinArm",
                         paste0("t_test_withinArm_", arm, "_B0_vs_", tp, ".csv"))
    write.csv(merged_ttest, file = outfile, row.names = FALSE)
    message("  T.test withinArm results saved: ", outfile)
  }
}