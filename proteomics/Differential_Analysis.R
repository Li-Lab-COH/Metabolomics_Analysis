# =============================================================================
# Proteomics Differential Expression Pipeline (Limma)
#
# - Load R libraries and import proteomics intensities, sample metadata, and annotations
# - Log2‐transform intensities and filter out the bottom 15% low‐variance proteins
# - Define helper functions for merging results with metadata and ensuring output directories
# - Intra‐arm (between diets) comparisons at D1–F12: unpaired IF vs NIF contrasts
# - Within‐arm (time‐course) analyses per diet: paired B0 vs D1–F12 with Patient_ID blocking
# - Compute and FDR‐correct statistics, retain proteins with raw P ≤ 0.05, merge with annotations
# - Save annotated CSV results into organized “intraArm/” and “withinArm/” subfolders
#
# Usage: source this script to execute all preprocessing, modeling, and output steps
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

#======================= Processing data =====================================
protData <- read.csv("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/ResultsTables_Filtered_1.csv")
metaData <- read_excel("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/meta_data_proteomics 2.xlsx")
prot_metadata <- read.csv("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/protein_metadata.csv")

# Preparing names
rownames(protData) <- protData$SeqID
protData$SeqID <- NULL

# Processing data
protData_log2 <- log2(protData)

# # Exploring data
# summary(unlist(protData[1032, ]))
# hist(unlist(protData[1032, ]), breaks = 20)
# hist(log2(unlist(protData[032, ])), breaks = 20)


#------------------- Removing low variance ---------------------------------

# # compute variances
# vars <- rowVars(as.matrix(protData_log2), na.rm=TRUE)
# 
# # pick a few candidate cutoffs
# cuts <- quantile(vars, probs = c(0.01, 0.05, 0.10, 0.15, 0.25))
# 
# # Base‐R histogram + vertical lines
# hist(vars,
#      breaks = 5000,
#      main   = "Protein variance distribution",
#      xlab   = "Variance",
#      xlim = c(0,0.3),
#      col    = "lightgray")
# abline(v = cuts, col = c("blue","green","orange","red"), lwd = 2)
# legend("topright",
#        legend = paste0(names(cuts)," = ", round(cuts,3)),
#        col    = c("blue","green","orange","red"),
#        lwd    = 2,
#        cex    = 0.8)

# Compute per‐protein variances
vars    <- rowVars(as.matrix(protData_log2), na.rm=TRUE)

# Establish the 15th-percentile cutoff
var_cut <- quantile(vars, 0.15)

# Subset your matrix to keep only proteins above that cutoff
dim(protData_log2)
protData_log2 <- protData_log2[vars > var_cut, ]
dim(protData_log2)
# Check dimensions
# cat("Before:", nrow(protData_log2), "proteins\n")
# cat("After: ", nrow(prot_clean),     "proteins (bottom 10% removed)\n")
#------------------------- Functions ----------------------------------

append_sig_results_with_metadata <- function(results_df, metadata_df, metadata_id_col = "SeqID") {
  
  # Process the results
  results_df <- results_df %>%
    filter(P.Value <= 0.05)
  results_df$ID <- rownames(results_df)
  
  # Merge metadata
  combined_df <- merge(results_df, metadata_df, by.x = "ID", by.y = metadata_id_col)
  
  # Rearrange columns: desired first, rest after (without duplication)
  desired_order <- c("Gene", "P.Value", "adj.P.Val", "logFC")
  remaining_cols <- setdiff(colnames(combined_df), desired_order)
  combined_df <- combined_df[, c(desired_order, remaining_cols)]
  
  return(combined_df)
}

# # Testing function
# sig_result <- append_sig_results_with_metadata(result, prot_metadata)
# colnames(sig_result)

ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}


#--------------------- Preparing folders --------------------------------------


# Define all output directories
dirs <- c(
  "~/1Work/RoseLab/Metabolomics/results_prot/difference/limma/",
  "~/1Work/RoseLab/Metabolomics/results_prot/difference/limma/withinArm/",
  "~/1Work/RoseLab/Metabolomics/results_prot/difference/limma/intraArm/"
)
sapply(dirs, ensure_dir)



#------------------------- loop ----------------------------------


# ================================
# LIMMA Analysis
# ================================

# IntraArm Comparisons (Between groups at the same time point)
# Make sure “Arm” is a factor with NIF first
metaData$Arm <- factor(metaData$Arm, levels = c("NIF", "IF"))

# Define exactly the time points you want to compare (exclude “B0” if you don’t need it here)
timepoints_of_interest <- c("D1", "D5", "F3", "F6", "F9", "F12")

for(tp in timepoints_of_interest) {
  message("Running limma intraArm analysis for time point: ", tp)
  
  # subset the metadata and the log2 data
  meta_tp <- metaData %>% filter(Time == tp)
  mat_tp  <- protData_log2[, meta_tp$Sample_ID, drop = FALSE]
  
  # check that both arms are present
  print(table(meta_tp$Arm))
  #    IF NIF 
  # e.g.  6   6
  
  # build design matrix
  design_tp <- model.matrix(~ Arm, data = meta_tp)
  colnames(design_tp) <- c("Intercept", "ArmIF")
  
  # inspect the first few rows of the design to make sure coding is correct
  print(head(cbind(meta_tp$Arm, design_tp)))
  #   meta_tp$Arm Intercept ArmIF
  # 1         NIF         1     0
  # 2         NIF         1     0
  # 3          IF         1     1
  # …         …          …     …
  
  # fit and extract the IF vs NIF contrast
  fit_tp     <- lmFit(mat_tp, design_tp) %>% eBayes()
  results_tp <- topTable(
    fit_tp,
    coef        = "ArmIF",
    adjust.method = "fdr",
    number      = Inf,
    sort.by     = "none"
  )
  
  # append your metadata (filtering for P ≤ 0.05 and re-ordering columns)
  processed_tp <- append_sig_results_with_metadata(
    results_df   = results_tp,
    metadata_df  = prot_metadata,
    metadata_id_col = "SeqID"
  )
  
  # write out the *processed* table
  outfile <- file.path(
    "~/1Work/RoseLab/Metabolomics/results_prot/difference/limma/intraArm/",
    paste0("limma_intraArm_", tp, "_IF_vs_NIF.csv")
  )
  write.csv(processed_tp, file = outfile, row.names = FALSE)
  message("Results saved: ", outfile)
}


#-------------------------Withinarm----------------------------------------

## (B) WithinArm Comparisons (Within each arm: later time points vs baseline B0)

# Create a Patient_ID column by pulling everything before the first underscore
metaData <- metaData %>%
  mutate(Patient_ID = str_extract(Sample_ID, "^[^_]+"))

# Confirm it worked
table(metaData$Patient_ID)
dim(metaData)
# IF01 IF02 IF05 ... 
#   6    5    6  ...

out_dir <- "~/1Work/RoseLab/Metabolomics/results_prot/difference/limma/withinArm/"

# Then in your within-Arm loop, block on that Patient_ID:


for(arm in c("NIF","IF")) {
  meta_arm <- filter(metaData, Arm == arm)
  
  for(tp in c("D1","D5","F3","F6","F9","F12")) {
    meta_arm_tp <- filter(meta_arm, Time %in% c("B0", tp))
    mat_tp      <- protData_log2[, meta_arm_tp$Sample_ID]
    
    meta_arm_tp$Time <- factor(meta_arm_tp$Time, levels = c("B0", tp))
    design <- model.matrix(~ Time, data = meta_arm_tp)
    colnames(design) <- c("Intercept", paste0("Time", tp))
    
    # Paired fit: block on Patient_ID
    corfit <- duplicateCorrelation(mat_tp,
                                   design,
                                   block = meta_arm_tp$Patient_ID)
    fit    <- lmFit(mat_tp,
                    design,
                    block       = meta_arm_tp$Patient_ID,
                    correlation = corfit$consensus)
    fit    <- eBayes(fit)
    
    res_paired <- topTable(fit,
                           coef          = paste0("Time", tp),
                           adjust.method = "fdr",
                           number        = Inf,
                           sort.by       = "none")
    
    sig_prots <- append_sig_results_with_metadata(res_paired,
                                                  prot_metadata)
    
    # Save
    outFile <- file.path(out_dir,
                         paste0("limma_withinArm_",
                                arm, "_B0_vs_", tp, ".csv"))
    write.csv(sig_prots, outFile, row.names = FALSE)
    message("    wrote: ", outFile)
   
  }
}






# 
# 
# #-------------------------Any values between 0 and 1? ------------------------
# 
# low_val <- protData[apply(protData, 1, function(row) any(row > 0 & row < 1)), ]
# # No, you can use log transformation safely
# 
# #------------------------- Testing loop -------------------------------------
# 

# 
# # Define timepoints of interest
# timepoints <- c("B0", "D1", "D5", "F3", "F6", "F9", "F12")
# arms <- unique(metaData$Arm)  # e.g., "IF", "NIF"
# 
# # Setting parameters to test loop
# a <- arms[1] #NIF
# tp <- timepoints[2] #D1
# 
# md_arm <- metaData %>% filter(Arm == a)
# 
# 
# ##### In time loop 
# comparison_samples <- md_arm %>%
#   filter(Time %in% c("B0", tp)) %>% 
#   pull(Sample_ID)
# 
# data_subset <- protData_log2[, comparison_samples, drop = FALSE]
# 
# 
# # Setting up column labels
# comparison_samples
# md_arm$Sample_ID
# 
# matches <- match(comparison_samples, md_arm$Sample_ID)
# matches
# md_arm$Time[matches]
# 
# time_group <- factor(md_arm$Time[matches],
#                      levels = timepoints)
# colnames(data_subset)
# metaData
# 
# #Setting up differential analysis
# # For each metabolite, fit a linear model where the expression depends on 
# # which timepoint the sample was collected at.
# 
# design <- model.matrix(~ 0 + time_group)
# colnames(design) <- levels(time_group)
# 
# design
# 
# # Build design matrix
# design <- model.matrix(~ 0 + time_group)
# colnames(design) <- levels(time_group)
# 
# # Contrast: current timepoint vs B0
# contrast.matrix <- makeContrasts(contrasts = paste(tp, "- B0"), levels = design)
# 
# # Fit model and compute statistics
# fit <- lmFit(data_subset, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# 
# # Get topTable with all metabolites
# result <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "none")
# 
# sig_result <- append_sig_results_with_metadata(result, prot_metadata)
# 
# # Save result to CSV
# file_name <- paste0("Differential_", a, "_", tp, "_vs_B0.csv")
# write.csv(result, file = file_name, row.names = TRUE)
# 
# # Progress message
# cat("Saved differential analysis for arm", a, "timepoint", tp, "vs B0 to", file_name, "\n")
# 



