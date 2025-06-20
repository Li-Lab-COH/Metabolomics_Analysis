# Load required libraries
library(limma)
library(dplyr)
BiocManager::install("limma")
# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12


metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID" )
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")

summary(unlist(metlData[9032, ]))
hist(unlist(metlData[9032, ]), breaks = 20)
hist(log2(unlist(metlData[9032, ])), breaks = 20)



#------------------------- Testing loop -------------------------------------

metlData_log2 <- log2(metlData + 1)

# Define timepoints of interest
timepoints <- c("B0", "D1", "D5", "F3", "F6", "F9", "F12")
arms <- unique(metaData$Arm)  # e.g., "IF", "NIF"

# Setting parameters to test loop
a <- arms[1] #NIF
tp <- timepoints[2] #D1

md_arm <- metaData %>% filter(Arm == a)


##### In time loop 
comparison_samples <- md_arm %>%
  filter(Time %in% c("B0", tp)) %>% 
  pull(Sample_ID)

data_subset <- metlData_log2[, comparison_samples, drop = FALSE]


# SEtting up column labels
comparison_samples
md_arm$Sample_ID

matches <- match(comparison_samples, md_arm$Sample_ID)
matches
md_arm$Time[matches]

time_group <- factor(md_arm$Time[matches],
        levels = timepoints)
colnames(data_subset)
metaData

#Setting up differential analysis
# For each metabolite, fit a linear model where the expression depends on which timepoint the sample was collected at.
design <- model.matrix(~ 0 + time_group)
colnames(design) <- levels(time_group)

design
# 
#   for (tp in setdiff(timepoints, "B0")) {
#     # Identify samples for this timepoint and baseline
#     comparison_samples <- md_arm %>%
#       filter(Time %in% c("B0", tp)) %>%
#       pull(Sample_ID)
# 
#     # Subset log2-transformed data for these samples
#     data_subset <- metlData_log2[, comparison_samples, drop = FALSE]
# 
#     # Make sure sample order in metadata matches column order in data
#     time_group <- factor(md_arm$Time[match(comparison_samples, md_arm$Sample_ID)],
#                          levels = timepoints)
    
    # Build design matrix
    design <- model.matrix(~ 0 + time_group)
    colnames(design) <- levels(time_group)
    
    # Contrast: current timepoint vs B0
    contrast.matrix <- makeContrasts(contrasts = paste(tp, "- B0"), levels = design)
    
    # Fit model and compute statistics
    fit <- lmFit(data_subset, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Get topTable with all metabolites
    result <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "none")
    
    # Save result to CSV
    file_name <- paste0("Differential_", a, "_", tp, "_vs_B0.csv")
    write.csv(result, file = file_name, row.names = TRUE)
    
    # Progress message
    cat("Saved differential analysis for arm", a, "timepoint", tp, "vs B0 to", file_name, "\n")
  }
}


#------------------------- suggested loop ----------------------------------



# -------------------------------
# Step 1: Prepare log2-transformed data
# -------------------------------

# Safely apply log2 transformation with a small offset to avoid log(0)
# Only apply if not already log-transformed
# You can adjust the pseudo count (1) if your values are very small

metlData_log2 <- log2(metlData + 1)

# Optional: Inspect distribution before and after log2 to verify
# hist(unlist(metlData[2000, ]), main = "Original")
# hist(unlist(metlData_log2[2000, ]), main = "Log2-transformed")

# -------------------------------
# Step 2: Differential analysis
# -------------------------------

# Define timepoints of interest
timepoints <- c("B0", "D1", "D5", "F3", "F6", "F9", "F12")
arms <- unique(metaData$Arm)  # e.g., "IF", "NIF"

# Loop over arms and timepoints to compare against baseline
for (a in arms) {
  # Subset metadata for current arm
  md_arm <- metaData %>% filter(Arm == a)
  
  for (tp in setdiff(timepoints, "B0")) {
    # Identify samples for this timepoint and baseline
    comparison_samples <- md_arm %>%
      filter(Time %in% c("B0", tp)) %>%
      pull(Sample_ID)
    
    # Subset log2-transformed data for these samples
    data_subset <- metlData_log2[, comparison_samples, drop = FALSE]
    
    # Make sure sample order in metadata matches column order in data
    time_group <- factor(md_arm$Time[match(comparison_samples, md_arm$Sample_ID)],
                         levels = timepoints)
    
    # Build design matrix
    #For each metabolite, fit a linear model where the expression depends on which timepoint the sample was collected at.
    design <- model.matrix(~ 0 + time_group)
    colnames(design) <- levels(time_group)
    
    # Contrast: current timepoint vs B0
    contrast.matrix <- makeContrasts(contrasts = paste(tp, "- B0"), levels = design)
    
    # Fit model and compute statistics
    fit <- lmFit(data_subset, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Get topTable with all metabolites
    result <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "none")
    
    # Save result to CSV
    file_name <- paste0("Differential_", a, "_", tp, "_vs_B0.csv")
    write.csv(result, file = file_name, row.names = TRUE)
    
    # Progress message
    cat("Saved differential analysis for arm", a, "timepoint", tp, "vs B0 to", file_name, "\n")
  }
}



# Testing QC values

mean(c(1117.08,	970.22,	1099.54,	912.8,	1240.92,	1001.48, 1033.5, 1018.93))
mean(c(5315.86, 5083.45, 5323.8, 5136.18, 4954.34, 4676.95, 5562, 5946.26))


# Z?
values = c(9217.31, 2962.35, 1479.78, 6062.47, 9947.16, 11192.67, 6397.76, 19920.43, 32211.57, 1887.68, 2322.79, 978.56, 2070.48, 6506.64, 2340.78, 13354.2, 2490.48, 7216.72, 5705.6, 756.79, 11268.81, 8882.67, 1588.32, 5333.98, 988.37, 4420.43, 5532.9, 6214.73, 18580.74, 19003.74, 8495.41, 7909.68, 16195.85, 16157.72, 14633.83, 7452.4, 7366.74, 6573.03, 7467.79, 7228.61, 7702.6, 9496.95, 10825.42, 8151.32, 16201.1, 12516.27, 4398.22, 5596.73, 4515.31, 3162.98, 4989.29, 4918.52, 833.96, 3840.72, 1246.28, 5722.09, 4157.88, 5328.78, 3454.98, 8832.4, 6254.04, 3845.41, 4456.84, 1282.77, 3970.84, 4431.15, 576.26, 9730.2, 838.03, 2507.56, 4334.64, 3581.32, 2361.17, 4804.25, 5713.96, 4852.78)
mean1 = mean(c(9217.31, 2962.35, 1479.78, 6062.47, 9947.16, 11192.67, 6397.76, 19920.43, 32211.57, 1887.68, 2322.79, 978.56, 2070.48, 6506.64, 2340.78, 13354.2, 2490.48, 7216.72, 5705.6, 756.79, 11268.81, 8882.67, 1588.32, 5333.98, 988.37, 4420.43, 5532.9, 6214.73, 18580.74, 19003.74, 8495.41, 7909.68, 16195.85, 16157.72, 14633.83, 7452.4, 7366.74, 6573.03, 7467.79, 7228.61, 7702.6, 9496.95, 10825.42, 8151.32, 16201.1, 12516.27, 4398.22, 5596.73, 4515.31, 3162.98, 4989.29, 4918.52, 833.96, 3840.72, 1246.28, 5722.09, 4157.88, 5328.78, 3454.98, 8832.4, 6254.04, 3845.41, 4456.84, 1282.77, 3970.84, 4431.15, 576.26, 9730.2, 838.03, 2507.56, 4334.64, 3581.32, 2361.17, 4804.25, 5713.96, 4852.78))
sd1 = sd(c(9217.31, 2962.35, 1479.78, 6062.47, 9947.16, 11192.67, 6397.76, 19920.43, 32211.57, 1887.68, 2322.79, 978.56, 2070.48, 6506.64, 2340.78, 13354.2, 2490.48, 7216.72, 5705.6, 756.79, 11268.81, 8882.67, 1588.32, 5333.98, 988.37, 4420.43, 5532.9, 6214.73, 18580.74, 19003.74, 8495.41, 7909.68, 16195.85, 16157.72, 14633.83, 7452.4, 7366.74, 6573.03, 7467.79, 7228.61, 7702.6, 9496.95, 10825.42, 8151.32, 16201.1, 12516.27, 4398.22, 5596.73, 4515.31, 3162.98, 4989.29, 4918.52, 833.96, 3840.72, 1246.28, 5722.09, 4157.88, 5328.78, 3454.98, 8832.4, 6254.04, 3845.41, 4456.84, 1282.77, 3970.84, 4431.15, 576.26, 9730.2, 838.03, 2507.56, 4334.64, 3581.32, 2361.17, 4804.25, 5713.96, 4852.78))


