# Load necessary libraries
library(limma)
library(dplyr)
# Timepoint definitions
# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12


metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID" )
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")
metlData_log2 <- log2(metlData + 1)

# Guiding question
# What happens over time in NIF (control)?
# 
# Does IF respond differently at each timepoint compared to NIF?

# -----------------------------
# 1. Data Subsetting & Preparation
# -----------------------------
# Define the timepoints to focus on
timepoints_of_interest <- c("B0", "D5", "F3", "F6")

# Subset the metadata to include only the desired timepoints and arms of interest
metaSubset <- metaData %>%
  filter(Time %in% timepoints_of_interest, Arm %in% c("IF", "NIF")) %>% 
  mutate(Patient_ID = gsub("_.*", "", Sample_ID))

# Ensure the samples in the metabolite data match the order in metaSubset
metlData_subset <- metlData_log2[, metaSubset$Sample_ID, drop = FALSE]

# Convert the Time variable to a factor with levels ordered as specified
# Levels defines the order of sequences, so the order is important when setting timepoints_of_interest
metaSubset$Time <- factor(metaSubset$Time, levels = timepoints_of_interest)



# Ensure Arm is a factor, using 'NIF' as the reference level
metaSubset$Arm <- factor(metaSubset$Arm, levels = c("NIF", "IF"))

# -----------------------------
# 2. Building the Design Matrix
# -----------------------------
# Create a design matrix with an interaction between Arm and Time.
# This will yield the following columns (by default):
#  - (Intercept) : expression for the reference group (Arm=NIF, Time=B0)
#  - ArmIF       : difference between arms at baseline (IF vs NIF)
#  - TimeD5, TimeF3, TimeF6 : differences from B0 for arm NIF
#  - ArmIF:TimeD5, ArmIF:TimeF3, ArmIF:TimeF6 : additional differences in the time effect for arm IF
design <- model.matrix(~ Arm * Time, data = metaSubset)
colnames(design) <- make.names(colnames(design))
print("Design matrix columns:")
print(colnames(design))

# -----------------------------
# 3. Fitting the Model with limma (Repeated Measures)
# -----------------------------
# If you have repeated measurements from the same patient, account for within-patient correlation.

corfit <- duplicateCorrelation(metlData_subset, design, block = metaSubset$Patient_ID)
fit <- lmFit(metlData_subset, design, block = metaSubset$Patient_ID, correlation = corfit$consensus)
fit <- eBayes(fit)

# -----------------------------
# 4. Defining Contrasts
# -----------------------------
# Here we define contrasts to extract:
# (a) Within-arm changes relative to baseline:
#     - For Arm NIF (baseline group): the contrast is given by the Time effect only.
#     - For Arm IF: the contrast is given by the Time effect plus the interaction (ArmIF:TimeX).
# (b) Between-arm differences at each timepoint.
# (c) The interaction terms (i.e. whether the time change differs between arms).
#
# Note: The column names used below are based on the model.matrix output.
#       They should be: "ArmIF", "TimeD5", "TimeF3", "TimeF6",
#                     "ArmIF:TimeD5", "ArmIF:TimeF3", etc.

contr.matrix <- makeContrasts(
  # Within-arm changes versus B0:
  NIF_D5_vs_B0 = TimeD5,                     # For NIF (control): change at D5 vs B0
  NIF_F3_vs_B0 = TimeF3,                     # For NIF: change at F3 vs B0
  NIF_F6_vs_B0 = TimeF6,                     # For NIF: change at F6 vs B0
  
  # For IF (treatment), the change is baseline + interaction:
  IF_D5_vs_B0 = TimeD5 + ArmIF.TimeD5,
  IF_F3_vs_B0 = TimeF3 + ArmIF.TimeF3,
  IF_F6_vs_B0 = TimeF6 + ArmIF.TimeF6,
  
  # Between-arm differences at each timepoint:
  IF_vs_NIF_B0 = ArmIF,                      # At baseline
  IF_vs_NIF_D5 = ArmIF + ArmIF.TimeD5,       # At D5
  IF_vs_NIF_F3 = ArmIF + ArmIF.TimeF3,       # At F3
  IF_vs_NIF_F6 = ArmIF + ArmIF.TimeF6,       # At F6
  
  # Interaction terms (difference in change over time between arms):
  Interact_D5 = ArmIF.TimeD5,
  Interact_F3 = ArmIF.TimeF3,
  Interact_F6 = ArmIF.TimeF6,
  
  levels = design
)

print("Contrast matrix:")
print(contr.matrix)

# Apply the contrasts to the fitted model
fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

# -----------------------------
# 5. Extracting and Saving Results
# -----------------------------
# Example: Extract results for the interaction at D5 (i.e. difference in D5 change between arms)
results_Interact_D5 <- topTable(fit2, coef = "Interact_D5", number = Inf)
print("Top results for interaction at D5 (difference in change between arms):")
print(results_Interact_D5)

# Similarly, extract results for other contrasts as needed:
results_IF_F3 <- topTable(fit2, coef = "IF_F3_vs_B0", number = Inf)
results_NIF_F3 <- topTable(fit2, coef = "NIF_F3_vs_B0", number = Inf)
results_between_F3 <- topTable(fit2, coef = "NIF_vs_IF_F3", number = Inf)

# Optionally, save some or all results to CSV files
write.csv(results_Interact_D5, file = "Interact_D5_results.csv", row.names = TRUE)
write.csv(results_IF_F3, file = "IF_F3_vs_B0_results.csv", row.names = TRUE)
write.csv(results_NIF_F3, file = "NIF_F3_vs_B0_results.csv", row.names = TRUE)
write.csv(results_between_F3, file = "Between_arms_F3_results.csv", row.names = TRUE)
