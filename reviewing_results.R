###################################
# 
# Author: Jonathan Anzules
# Purpose: Analyzing results from the limma time series analysis performed in "Limma_time_series_analysis.R"
# Input: Limma results and metabolite cleaned dataset
# Output: Summary figures related to the results, manually saved, so becareful of overwritting


# Loading Packages
packages <- c("dplyr", "tidyr", "ggplot2", "rlang", "limma")
lapply(packages, library, character.only = TRUE)



#----------------------- Preparing input data -----------------------------

fit_results <- readRDS("~/Roselab/Metabolite/results/time_series_analysis/fit2_results_B0toF6_filtered.rds")
metabolite_metadata <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_metadata.csv")

metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")


# Preparing metabolite data
metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID" )
metlData_log2 <- log2(metlData + 1)
metlData_log2 <- metlData_log2 %>% 
  mutate(
    Index = sub("\\|.*", "", rownames(metlData_log2))
  )

#-------------------- Extract Function n Dir -----------------------------------
extractContrastResults <- function(fit, coef_name, pValue = 0.05, logfc_cutoff = 0, 
                                   annotate = NULL, id_col = "Index",
                                   significant_only = FALSE, incl_metl_data = FALSE,
                                   metl_data = NULL) {
  # Load needed lib
  require(dplyr)
  
  # Extract results from fit2
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  
  # Add significance and direction labels
  res <- res %>%
    mutate(
      Index = sub("\\|.*", "", rownames(res)),
      Significant = P.Value < pValue,
      Direction = case_when(
        logFC > logfc_cutoff ~ "Up",
        logFC < -logfc_cutoff ~ "Down",
        TRUE ~ "No change"
      ),
      Contrast = coef_name,
    )
  
  if(significant_only){
    res <- res %>% filter(Significant)
  }
  
  # Merge in metabolite annotations, if provided
  if (!is.null(annotate)) {
    res <- left_join(res, annotate, by = id_col)
  }
  
  # Merge metabolite data
  if(incl_metl_data){
    res <- left_join(res, metl_data, by = id_col)
  }
  
  rownames(res) <- res$Index
  res$Significant <- NULL
  return(res)
}

# Function to create direcotries
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}

#----------------------- Extracting Data --------------------------------------

# Checklist of plots
# Interact_D5 - DONE
# Interact_F3 - DONE

test_extract <- extractContrastResults(fit_results, "Interact_D5",
                                       annotate = metabolite_metadata,
                                       significant_only = TRUE,
                                       incl_metl_data = TRUE,
                                       metl_data = metlData_log2
                                      )

# test_extract <- extractContrastResults(fit_results, "Interact_F3",
#                                        annotate = metabolite_metadata,
#                                        significant_only = TRUE,
#                                        incl_metl_data = TRUE,
#                                        metl_data = metlData_log2
# )


#------------ Preparing directories and variables for plotting-----------------
# timepoints_of_interest <- c("B0", "D5", "F3", "F6")
timepoints_of_interest <- c("B0", "D5", "F3")

##### What lf change are we looking at?
test_extract <- subset(test_extract, logFC > 0)
# test_extract <- subset(test_extract, logFC < 0)

##### What class are we looking at
# group_by_class <- sym("Class.I") # Class.I or Class.II
group_by_class <- sym("Class.II") # Class.I or Class.II

# Preparing directories
# TODO: add an if check depending on whether subsetted for positive or negative lfchange
out_dir <- "~/Roselab/Metabolite/results/time_series_analysis/Interact_D5/posLF/ClassII"
ensure_dir(out_dir)

#For combined summaries
out_sum <- "~/Roselab/Metabolite/results/time_series_analysis/Interact_D5/LFSummary/ClassII"
ensure_dir(out_sum)






#------------------------ Plotting trends -------------------------------------

# List the known metadata columns (which do NOT contain sample expression data).
metadata_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
                   "Index", "Direction", "Contrast",
                   "Compounds", "Class.I", "Class.II", "Formula", "Mode",
                   "Molecular.weight..Da.", "RT..min.", "Adduct", "Mass.error",
                   "Level", "score", "LC.mode", "Endogenous.Exogenous",
                   "CAS", "PubChem.CID", "HMDB", "Metlin", "cpd_ID", 
                   "QC01", "QC02", "QC03", "QC04", "QC05", "QC06", "QC07", "QC08",
                   "kegg_map")

# The remaining columns will be sample IDs.
sample_cols <- setdiff(colnames(test_extract), metadata_cols)

# --- Pivot the data into long format ---
# Each row will correspond to one metabolite in one sample.
test_extract_long <- test_extract %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample_ID",
    values_to = "Expression"
  )

#### Join with the metadata
# 'metaData' contains information like Sample_ID, Time, Arm, etc.
#  metaData$Sample_ID uniquely identifies each sample and matches those in test_extract_long.
combined_data <- test_extract_long %>%
  inner_join(metaData, by = "Sample_ID")

#Filter by timepoints of interest
combined_data <- combined_data %>% 
  filter(Time %in% timepoints_of_interest)

#### Group by metabolite class and time, then summarize

trend_summary <- combined_data %>%
  group_by(!!group_by_class, Time, Arm) %>% 
  summarize(
    mean_expr = mean(Expression, na.rm = TRUE),
    sd_expr   = sd(Expression, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  )

# --------------------- Begin plotting -------------------------------------
# Separate facets for each class, with lines for each Arm (IF and NIF)
p1 <- ggplot(trend_summary, aes(x = Time, y = mean_expr, group = Arm, color = Arm)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
                width = 0.2) +
  facet_wrap(vars(!!group_by_class), scales = "free_y") +
  labs(title = "Metabolite Class Trends Over Time",
       x = "Timepoint",
       y = "Mean log2 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
# print(p1)



p1.2 <- ggplot() +
  # Plot individual sample points
  geom_point(data = combined_data, 
             aes(x = Time, y = Expression, color = Arm),
             position = position_jitter(width = 0.1), alpha = 0.5, size = 0.5) +
  
  # Plot mean trends with error bars
  geom_line(data = trend_summary, 
            aes(x = Time, y = mean_expr, group = Arm, color = Arm), 
            size = 1.2) +
  
  geom_point(data = trend_summary, 
             aes(x = Time, y = mean_expr, color = Arm), 
             size = 2.5) +
  
  geom_errorbar(data = trend_summary, 
                aes(x = Time, y = mean_expr, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr, color = Arm),
                width = 0.15) +
  
  # Facet by the selected class
  facet_wrap(vars(!!group_by_class), scales = "free_y") +
  
  labs(title = "Metabolite Class Trends Over Time with Individual Points",
       x = "Timepoint",
       y = "log2 Expression (individual + mean ± SD)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# print(p1.2)


# see separate plots for the two Arms for a given class, facet by both Class and Arm.
p2 <- ggplot(trend_summary, aes(x = Time, y = mean_expr, group = 1)) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 2, color = "blue") +
  geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
                width = 0.2) +
  facet_grid(rows = vars(Arm), cols = vars(!!group_by_class), scales = "free_y") +
  labs(title = "Metabolite Class Trends by Arm",
       x = "Timepoint",
       y = "Mean log2 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
print(p2)



# Violin Plots
p_violin <- ggplot(test_extract, aes(x = !!group_by_class, y = logFC, fill = !!group_by_class)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "black") +
  labs(title = "Violin Plot of Interact_D5 LogFC by Metabolite Class",
       subtitle = "Positive values: IF shows a larger change (B0 -> D5) than NIF",
       x = "Metabolite Class",
       y = "Log2 Fold Change (Interact_D5)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
print(p_violin)


#---------------------- Summarizing data by class ----------------------------

# Summarize the data by class using the grouping variable
class_summary <- test_extract %>%
  group_by(!!group_by_class) %>%
  summarize(
    mean_logFC = mean(logFC, na.rm = TRUE),
    se_logFC   = sd(logFC, na.rm = TRUE)/sqrt(n()),
    count      = n(),
    .groups    = "drop"
  )

# TODO: where there is a (), when I make a function for all of this, I should add
#       a description as to where the data is from, I should do this for every plot... ?
p_bar <- ggplot(class_summary, aes(x = !!group_by_class, y = mean_logFC, fill = !!group_by_class)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_logFC - se_logFC, ymax = mean_logFC + se_logFC), width = 0.2) +
  labs(title = "Mean () LogFC by Metabolite Class",
       subtitle = "Error bars represent standard error; note that mixed positive & negative values may cancel",
       x = "Metabolite Class",
       y = "Mean Log2 Fold Change ()") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
# print(p_bar)

# bile_acids <- subset(test_extract, Class.I == "Bile acids")

p_box <- ggplot(test_extract, aes(x = !!group_by_class, y = logFC, fill = !!group_by_class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # outlier.shape = NA suppresses default outliers
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "black") +
  labs(title = "Boxplot of Interact_D5 LogFC by Metabolite Class",
       subtitle = "Boxplot shows median and IQR; deviations indicate heterogeneity",
       x = "Metabolite Class",
       y = "Log2 Fold Change (Interact_D5)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
# print(p_box)


################################################################################
# ONLY DO THIS IF YOU HAVE NOT SEPERATED BY LF CHANGE
################################################################################
# Summarize counts by class and direction
direction_summary <- test_extract %>%
  group_by(!!group_by_class, Direction) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n))

p_stacked <- ggplot(direction_summary, aes(x = !!group_by_class, y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Counts of Up/Down/No Change Metabolites by Class",
       subtitle = "Direction based on logFC cutoff thresholds",
       x = "Metabolite Class",
       y = "Number of Metabolites") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
# print(p_stacked)


# Proportions
p_stacked_prop <- ggplot(direction_summary, aes(x = !!group_by_class, y = prop, fill = Direction)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(title = "Proportion of Up/Down/No Change Metabolites by Class",
       subtitle = "Staked to 100%; higher proportion indicates consistency",
       x = "Metabolite Class",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
# print(p_stacked_prop)

#------------------------- Saving all figures -------------------------------

ggsave(filename = file.path(out_dir, "class.Trends.png" ), plot = p1, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "class.Trends.points.png" ), plot = p1.2, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "class.Trends.separateArms.png" ), plot = p2, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "LFC.violin.png" ), plot = p_violin, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "LFC.Barchart.png" ), plot = p_bar, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "LFC.Boxplot.png" ), plot = p_box, width = 10, height = 6, dpi = 300)


#LOG FOLD CHANGE SUMMARY PLOTS IF NO SUBSETTING
ggsave(filename = file.path(out_sum, "LFC.stackedBar.png" ), plot = p_stacked, width = 10, height = 6, dpi = 300)
ggsave(filename = file.path(out_sum, "LFC.stackedBar.Proportions.png" ), plot = p_stacked_prop, width = 10, height = 6, dpi = 300)



