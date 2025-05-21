#################################
# Load data
#################################
# Imputed dataset
metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID")
# Human metadata
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")
# Metabolite metadata
metabolite_meta_data <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_metadata.csv")
# Metabolites to remove
metabolites_to_remove <- scan(
  "~/Roselab/Metabolite/data/data_for_analysis/metabolites_to_remove.txt",
  what = character(),
  sep  = "\n",
  quiet= TRUE
)



#################################
# Function
#################################

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

#################################
# T testing
# TODO: Need to add the main testing into a larger loop to test all time points
#################################

intra_ttest_results <- list()
timepoints_of_interest_ttest <- c("D5", "F3", "F6")

meta_tp <- metaSubset[metaSubset$Time == "F6", ]
metlData_tp <- metlData_log2_filtered[, meta_tp$Sample_ID, drop = FALSE]
group_vec <- factor(as.character(meta_tp$Arm), levels = c("NIF", "IF"))

# Extract metabolite values (row 1) and make it numeric
met_values <- as.numeric(metlData_tp[1, ])
names(met_values) <- colnames(metlData_tp)  # keep sample IDs

# Split into two groups manually
group1 <- met_values[group_vec == "NIF"]
group2 <- met_values[group_vec == "IF"]

# Run the unpaired t-test
tt <- t.test(group1, group2, paired = FALSE)


# This is the prefered way to do a t.test, but because some metabolites cause a failure.,
# The loop helps identify which metabolites are causing the problems.
# pvals <- sapply(1:nrow(metlData_tp), function(i) {
#   t.test( metlData_tp[i, group_vec=="NIF"],
#           metlData_tp[i, group_vec=="IF"],
#           paired = FALSE )$p.value
# })


results <- data.frame(
  Metabolite = character(),
  t.stat     = numeric(),
  p.value    = numeric(),
  stringsAsFactors = FALSE
)

# loop over each metabolite (row)
for (i in seq_len(nrow(metlData_tp))) {
  met_name <- rownames(metlData_tp)[i]
  
  # Skip if in the removal list
  if (met_name %in% metabolites_to_remove) {
    cat("Skipping:", met_name, "(in removal list)\n")
    next
  }
  
  
  group1   <- metlData_tp[i, group_vec == "NIF"]
  group2   <- metlData_tp[i, group_vec == "IF"]
  
  # print the metabolite name
  cat("----- Testing:", met_name, "-----\n")
  
  # print the raw values
  cat("  NIF:", paste0(round(group1, 3), collapse = ", "), "\n")
  cat("   IF:", paste0(round(group2, 3), collapse = ", "), "\n")
  
  # run the t-test
  tt <- t.test(group1, group2, paired = FALSE)
  cat(sprintf("  t = %.3f, p = %.4g\n\n", tt$statistic, tt$p.value))
  
  # append to results
  results <- rbind(
    results,
    data.frame(
      Metabolite = met_name,
      t.stat     = unname(tt$statistic),
      p.value    = tt$p.value,
      stringsAsFactors = FALSE
    #TODO: Need to use the process_results function here to add the results back to metabolite metadata
    #TODO: results should be saved in appropriate folder, depending on loop iteration
    )
  )
}


# TODO: FDR‐corrected column - Should test different methods, but nothing seems significant in terms
# of adjusted p-value.
results$FDR <- p.adjust(results$p.value, method = "fdr")



#################################
# TODO: Need to add the a within arm testing (e.g. D5 vs F3; another loop section - see Statistical_testing.R )
#################################




