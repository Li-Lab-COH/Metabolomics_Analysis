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



# pvals <- sapply(1:nrow(metlData_tp), function(i) {
#   t.test( metlData_tp[i, group_vec=="NIF"],
#           metlData_tp[i, group_vec=="IF"],
#           paired = FALSE )$p.value
# })
metabolites_to_remove <- scan(
  "~/Roselab/Metabolite/data/data_for_analysis/metabolites_to_remove.txt",
  what = character(),
  sep  = "\n",
  quiet= TRUE
)

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
  
  # 1) print the metabolite name
  cat("----- Testing:", met_name, "-----\n")
  # 2) print the raw values
  cat("  NIF:", paste0(round(group1, 3), collapse = ", "), "\n")
  cat("   IF:", paste0(round(group2, 3), collapse = ", "), "\n")
  
  # 3) run the t-test
  tt <- t.test(group1, group2, paired = FALSE)
  cat(sprintf("  t = %.3f, p = %.4g\n\n", tt$statistic, tt$p.value))
  
  # 4) append to results
  results <- rbind(
    results,
    data.frame(
      Metabolite = met_name,
      t.stat     = unname(tt$statistic),
      p.value    = tt$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# finally, add an FDR‐corrected column
results$FDR <- p.adjust(results$p.value, method = "fdr")

# take a peek
head(results)
