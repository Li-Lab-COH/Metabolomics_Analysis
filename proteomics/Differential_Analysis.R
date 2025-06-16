# Load required libraries
library(limma)
library(dplyr)
library(readxl)

# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12


protData <- read.csv("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/ResultsTables_Filtered_1.csv")
metaData <- read_excel("~/1Work/RoseLab/Metabolomics/data/proteomics/for_analysis/meta_data_proteomics 2.xlsx")

rownames(protData) <- protData$SeqID
protData$SeqID <- NULL
  


summary(unlist(protData[1032, ]))
hist(unlist(protData[1032, ]), breaks = 20)
hist(log2(unlist(protData[032, ])), breaks = 20)

#-------------------------Any values between 0 and 1? ------------------------

low_val <- protData[apply(protData, 1, function(row) any(row > 0 & row < 1)), ]
# No, you can use log transformation safely

#------------------------- Testing loop -------------------------------------

protData_log2 <- log2(protData)

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

data_subset <- protData_log2[, comparison_samples, drop = FALSE]


# Setting up column labels
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
#     data_subset <- protData_log2[, comparison_samples, drop = FALSE]
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
#   }
# }


#------------------------- suggested loop ----------------------------------



# -------------------------------
# Step 1: Prepare log2-transformed data
# -------------------------------

# Safely apply log2 transformation with a small offset to avoid log(0)
# Only apply if not already log-transformed
# You can adjust the pseudo count (1) if your values are very small

protData_log2 <- log2(protData + 1)

# Optional: Inspect distribution before and after log2 to verify
# hist(unlist(protData[2000, ]), main = "Original")
# hist(unlist(protData_log2[2000, ]), main = "Log2-transformed")

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
    data_subset <- protData_log2[, comparison_samples, drop = FALSE]
    
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
