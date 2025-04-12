metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID")
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")
# Log2 transform (adding 1 to avoid log(0))
metlData_log2 <- log2(metlData + 1)
colnames(metlData_log2)

library(ggplot2)
library(dplyr)

plot_metabolite_trends <- function(metabolite_names, output_dir, metlData_log2, metaData) {
  
  if (!dir.exists(output_dir)){
    print("Folder created! <3")
    dir.create(output_dir, recursive = TRUE)
  } else{
    print("Folder exists! <3")
  }
  
  # Ensure that metabolite_names is a character vector (even for single values)
  metabolite_names <- as.character(metabolite_names)
  
  # Define the timepoints of interest (using "B0" as baseline)
  valid_times <- c("B0", "D5", "F3", "F6")
  
  for (met in metabolite_names) {
    
    # Check that the metabolite exists in the data
    if (!met %in% rownames(metlData_log2)) {
      warning(paste("Metabolite", met, "not found in metlData_log2"))
      next
    }
    
    # Extract metabolite values; the names of the values are the sample IDs.
    metabolite_values <- metlData_log2[met, ]
    print(metabolite_values)
    df_met <- data.frame(Sample_ID = colnames(metlData_log2),
                         Value = as.numeric(metabolite_values),
                         stringsAsFactors = FALSE)
    
    # Merge with metadata on the common "Sample_ID" column.
    df_plot <- merge(df_met, metaData, by = "Sample_ID")
    
    # Filter to include only samples from the valid timepoints.
    df_plot <- df_plot %>% filter(Time %in% valid_times)
    if(nrow(df_plot) == 0){
      warning(paste("No data for metabolite", met, "for the specified timepoints"))
      next
    }
    
    # Use the original Time column directly.
    # Setting factor levels ensures that the x-axis appears in the desired order.
    df_plot$Time <- factor(df_plot$Time, levels = valid_times)
    
    # Compute summary statistics per Time and Arm group.
    df_summary <- df_plot %>%
      group_by(Time, Arm) %>%
      summarise(meanValue = mean(Value, na.rm = TRUE),
                sdValue = sd(Value, na.rm = TRUE),
                n = n(),
                seValue = sdValue / sqrt(n),
                .groups = "drop")
    
    print(df_summary)
    
    # Define a dodge width for grouping values by Arm.
    dodge_width <- 0.3
    p <- ggplot(df_plot, aes(x = Time, y = Value, color = Arm)) +
      geom_jitter(position = position_dodge(width = dodge_width), alpha = 0.7, size = 2) +
      geom_point(data = df_summary,
                 aes(x = Time, y = meanValue, color = Arm, group = Arm),
                 size = 4, shape = 18, position = position_dodge(width = dodge_width),
                 inherit.aes = FALSE) +
      geom_errorbar(data = df_summary,
                    aes(x = Time, ymin = meanValue - seValue, ymax = meanValue + seValue, color = Arm, group = Arm),
                    width = 0.2, position = position_dodge(width = dodge_width),
                    inherit.aes = FALSE) +
      labs(title = met,
           x = "Time (B0, D5, F3, F6)",
           y = "Log₂-transformed abundance") +
      theme_bw()
    
    # Define the filename based on the metabolite name.
    filename <- file.path(output_dir, paste0(met, "dotPlot", ".png"))
    
    # Save the plot (adjust width, height, or file type as needed)
    ggsave(filename, plot = p, width = 8, height = 6)
    
    message("Plot saved for metabolite: ", met, " at ", filename)
  }
}



plot_metabolite_boxplots <- function(metabolite_names, output_dir, metlData_log2, metaData) {
  
  # Ensure metabolite_names is a character vector (for single or multiple metabolites)
  metabolite_names <- as.character(metabolite_names)
  
  # Define timepoints to include (with "B0" as baseline)
  valid_times <- c("B0", "D5", "F3", "F6")
  
  # Loop over each metabolite provided
  for (met in metabolite_names) {
    
    # Verify that the metabolite exists
    if (!met %in% rownames(metlData_log2)) {
      warning(paste("Metabolite", met, "not found in metlData_log2"))
      next
    }
    
    # Extract the metabolite values (the column names are your sample IDs)
    metabolite_values <- metlData_log2[met, ]
    df_met <- data.frame(Sample_ID = colnames(metlData_log2),
                         Value = as.numeric(metabolite_values),
                         stringsAsFactors = FALSE)
    
    # Merge these values with the metadata using the Sample_ID column
    df_plot <- merge(df_met, metaData, by = "Sample_ID")
    
    # Filter to include only the desired timepoints
    df_plot <- df_plot %>% filter(Time %in% valid_times)
    if(nrow(df_plot) == 0){
      warning(paste("No data for metabolite", met, "for the specified timepoints"))
      next
    }
    
    # Convert Time to a factor to ensure proper ordering on the x-axis
    df_plot$Time <- factor(df_plot$Time, levels = valid_times)
    
    # Define a dodge width for separating the boxplots by Arm at each Time point
    dodge_width <- 0.3
    
    # Create the boxplot. The 'fill' aesthetic will separate the two Arms
    p <- ggplot(df_plot, aes(x = Time, y = Value, fill = Arm)) +
      geom_boxplot(position = position_dodge(width = dodge_width)) +
      # Optionally overlay jittered points to visualize individual data points:
      geom_jitter(aes(color = Arm), 
                  position = position_dodge(width = dodge_width), 
                  alpha = 0.7, size = 2) +
      labs(title = met,
           x = "Time (B0, D5, F3, F6)",
           y = "Log₂-transformed abundance") +
      theme_bw()
    
    # Create a filename from the metabolite name
    filename <- file.path(output_dir, paste0(met, "boxplot", ".png"))
    
    # Save the plot (adjust dimensions or file type as needed)
    ggsave(filename, plot = p, width = 8, height = 6)
    
    message("Box-and-whiskers plot saved for metabolite: ", met, " at ", filename)
  }
}

# 
# output_directory <- "~/Roselab/Metabolite/results/Figures/MetabolitePlots/"
# 
# 
# plot_metabolite_boxplots(metabolite_names = c("MW0113847|Acetildenafil"),
#                          output_dir = output_directory,
#                          metlData_log2 = metlData_log2,
#                          metaData = metaData)
# 
# plot_metabolite_trends(metabolite_names = "MW0113847|Acetildenafil",
#                        output_dir = output_directory,
#                        metlData_log2 = metlData_log2,
#                        metaData = metaData)


#-------------------------- Common metabolite plots - peptides ----------------

# outDir_common_peptides <- "~/Roselab/Metabolite/results/Figures/MetabolitePlots/common_metabolites/peptides/"
# 
# #Pepties
# common_pepties <- c("MW0151607|Ile-Thr-Tyr-Asp", "MW0159127|Val-Phe-Phe-Asn-Gly",
#                     "MW0151164|His-HoPhe-OH", "MW0158431|Tyr-Gln-Asp")
# 
# plot_metabolite_trends(metabolite_names = common_pepties,
#                        output_dir = outDir_common_peptides,
#                        metlData_log2 = metlData_log2,
#                        metaData = metaData)
# 
# plot_metabolite_boxplots(metabolite_names = common_pepties,
#                          output_dir = outDir_common_peptides,
#                          metlData_log2 = metlData_log2,
#                          metaData = metaData)


# ------------------ Common metabolites - hormones -----------------------------

outDir_common_hormones <- "~/Roselab/Metabolite/results/Figures/MetabolitePlots/common_metabolites/hormones/"

common_hormone <- c("MW0017085|Cholesterol sulfate", "MW0126450|Riboflavin reduced",
                    "MW0063752|Ursocholyltaurine", "MW0147772|Coproporphyrin I",
                    "MW0105369|6-Sulfatoxymelatonin", "MW0010768|Cloprostenol"
                    )

plot_metabolite_trends(metabolite_names = common_hormone,
                       output_dir = outDir_common_hormones,
                       metlData_log2 = metlData_log2,
                       metaData = metaData)

plot_metabolite_boxplots(metabolite_names = common_hormone,
                         output_dir = outDir_common_hormones,
                         metlData_log2 = metlData_log2,
                         metaData = metaData)















