# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12

metlData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = "ID")
metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv")
prot_results <- read.csv("~/Roselab/Metabolite/results/results_all_fromComp/results_prot/intraArm/limma_intraArm_ALL_timepoints_IF_vs_NIF_combined.csv")

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



plot_metabolite_boxplots <- function(
    metabolite_names, 
    output_dir,
    valid_times = c("B0", "D1", "D5", "F3", "F6"),
    metlData_log2, 
    metaData) {
  
  # Ensure metabolite_names is a character vector (for single or multiple metabolites)
  metabolite_names <- as.character(metabolite_names)
  
  # Define timepoints to include (with "B0" as baseline)
  # valid_times <- c("B0", "D1", "D5", "F3", "F6")
  
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
    df_plot <- merge(df_met, metaData, by = "Sample_ID", all.x = TRUE) %>%
      filter(Time %in% valid_times) %>%
      filter(!grepl("^MB", Patient_ID))
    
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

library(ggplot2)
library(dplyr)
library(stringr)

plot_metabolite_patient_trends <- function(
    metabolite_names,
    output_dir,
    metlData_log2,
    metaData,
    valid_times = c("B0", "D1", "D5", "F3", "F6"),
    facet_by_arm = TRUE,
    arm_colors = c("IF"  = "#E15759",  # reddish
                   "NIF" = "#4E79A7"), # blue
    plot_title = NULL,
    verbose_duplicates = TRUE
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created folder: ", output_dir)
  } else {
    message("Folder exists: ", output_dir)
  }
  
  metabolite_names <- as.character(metabolite_names)
  
  # column checks
  need_cols <- c("Sample_ID","Time","Arm","Patient_ID")
  miss <- setdiff(need_cols, names(metaData))
  if (length(miss)) stop("metaData missing columns: ", paste(miss, collapse = ", "))
  
  metaData <- metaData %>%
    mutate(Arm = factor(Arm, levels = names(arm_colors)))
  
  for (met in metabolite_names) {
    if (!met %in% rownames(metlData_log2)) {
      warning(sprintf("Metabolite '%s' not found in metlData_log2; skipping.", met))
      next
    }
    
    vals <- metlData_log2[met, ]
    df_met <- data.frame(
      Sample_ID = colnames(metlData_log2),
      Value     = as.numeric(vals),
      stringsAsFactors = FALSE
    )
    
    # merge, keep valid times, drop MB patients
    df_plot_raw <- merge(df_met, metaData, by = "Sample_ID", all.x = TRUE) %>%
      filter(Time %in% valid_times) %>%
      filter(!grepl("^MB", Patient_ID))
    
    if (nrow(df_plot_raw) == 0) {
      warning(sprintf("No rows for '%s' after filtering (valid_times/MB); skipping.", met))
      next
    }
    
    # detect duplicates at Patient_ID x Time x Arm
    dup_tbl <- df_plot_raw %>%
      count(Patient_ID, Time, Arm, name = "n_obs") %>%
      filter(n_obs > 1)
    
    if (verbose_duplicates && nrow(dup_tbl) > 0) {
      message("Averaging duplicates for '", met, "' at these Patient/Time/Arm combos:")
      print(dup_tbl)
    }
    
    # average duplicates -> one row per Patient_ID x Time x Arm
    df_plot <- df_plot_raw %>%
      group_by(Patient_ID, Time, Arm) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        Time = factor(Time, levels = valid_times),
        Arm  = droplevels(Arm)
      )
    
    if (nrow(df_plot) == 0) {
      warning(sprintf("No rows for '%s' after aggregation; skipping.", met))
      next
    }
    
    # title
    this_title <- if (is.null(plot_title)) met else plot_title
    
    # plot: one line per patient, colored by Arm
    p <- ggplot(
      df_plot,
      aes(x = Time, y = Value, group = Patient_ID, color = Arm)
    ) +
      geom_line(alpha = 0.55, linewidth = 0.9) +
      geom_point(alpha = 0.9, size = 2) +
      scale_color_manual(values = arm_colors, drop = FALSE, name = "Cohort") +
      labs(
        title = this_title,
        x = "Time",
        y = expression("Log"[2]*" abundance")
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "top"
      )
    
    if (isTRUE(facet_by_arm)) {
      p <- p + facet_wrap(~ Arm, nrow = 1, scales = "fixed")
    }
    
    safe_name <- stringr::str_replace_all(met, "[^A-Za-z0-9._-]", "_")
    out_png   <- file.path(output_dir, paste0(safe_name, "_patientLines.png"))
    ggsave(out_png, plot = p, width = 9, height = 5.5, dpi = 300)
    message("Saved: ", out_png)
  }
}


plot_metabolite_patient_trends_overlap <- function(
    metabolite_names,
    output_dir,
    metlData_log2,
    metaData,
    valid_times = c("B0","D1", "D5", "F3", "F6"),
    arm_colors = c("IF"  = "#E15759",
                   "NIF" = "#4E79A7"),
    plot_title = NULL,
    verbose_duplicates = TRUE,
    point_jitter_width = 0.12   # set to 0 for no jitter
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  metabolite_names <- as.character(metabolite_names)
  need_cols <- c("Sample_ID","Time","Arm","Patient_ID")
  miss <- setdiff(need_cols, names(metaData))
  if (length(miss)) stop("metaData missing columns: ", paste(miss, collapse = ", "))
  
  metaData <- metaData %>% mutate(Arm = factor(Arm, levels = names(arm_colors)))
  
  for (met in metabolite_names) {
    if (!met %in% rownames(metlData_log2)) {
      warning(sprintf("Metabolite '%s' not found; skipping.", met)); next
    }
    
    vals <- metlData_log2[met, ]
    df_met <- data.frame(Sample_ID = colnames(metlData_log2),
                         Value = as.numeric(vals), stringsAsFactors = FALSE)
    
    # merge, keep valid times, drop MB patients
    df_plot_raw <- merge(df_met, metaData, by = "Sample_ID", all.x = TRUE) %>%
      filter(Time %in% valid_times) 
    # %>%
    #   filter(!grepl("^MB", Patient_ID))
    
    if (!nrow(df_plot_raw)) { warning(sprintf("No rows for '%s' after filtering.", met)); next }
    
    # detect duplicates at Patient×Time×Arm
    dup_tbl <- df_plot_raw %>% count(Patient_ID, Time, Arm, name = "n_obs") %>% filter(n_obs > 1)
    if (verbose_duplicates && nrow(dup_tbl) > 0) {
      message("Averaging duplicates for '", met, "':"); print(dup_tbl)
    }
    
    # average duplicates -> one row per Patient×Time×Arm (used for both lines and points)
    df_avg <- df_plot_raw %>%
      group_by(Patient_ID, Time, Arm) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      mutate(Time = factor(Time, levels = valid_times),
             Arm  = droplevels(Arm))
    
    if (!nrow(df_avg)) { warning(sprintf("No rows for '%s' after aggregation.", met)); next }
    
    this_title <- if (is.null(plot_title)) met else plot_title
    
    p <- ggplot() +
      # per-patient trajectories (averaged)
      geom_line(
        data = df_avg,
        aes(x = Time, y = Value, group = Patient_ID, color = Arm),
        alpha = 0.55, linewidth = 0.9
      ) +
      # single point per Patient×Time×Arm (same averaged value as the line)
      geom_point(
        data = df_avg,
        aes(x = Time, y = Value, color = Arm),
        position = position_jitter(width = point_jitter_width, height = 0),
        alpha = 0.9, size = 2
      ) +
      scale_color_manual(values = arm_colors, drop = FALSE, name = "Cohort") +
      labs(title = this_title, x = "Time", y = expression("Log"[2]*" abundance")) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "top") #+
      #expand_limits(y = 0)
    
    out_png <- file.path(output_dir, paste0(stringr::str_replace_all(met, "[^A-Za-z0-9._-]", "_"),
                                            "_patientLines.png"))
    ggsave(out_png, plot = p, width = 9, height = 5.5, dpi = 300)
    message("Saved: ", out_png)
  }
}

plot_heatmap <- function(
    protein_names = c("prot1","prot2"),
    output_dir = ".",
    valid_times = c("B0","D1","D5","F3","F6","F9","F12"),
    prot_results,
    zscore_rows = FALSE,
    fname = "heatmap_proteins.png"
){
  # Dependencies
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install pheatmap: install.packages('pheatmap')")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install stringr")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Standardize column names
  df <- prot_results
  # Try to map possible column names
  col_symbol <- c("symbol","Symbol","Gene","gene","GeneSymbol","Protein")
  col_time   <- c("timepoint","Timepoint","time","Time","Label","contrast")
  col_logfc  <- c("logFC","LogFC","log2FC","log2FoldChange","LFC")
  
  map_col <- function(df, candidates){
    nm <- names(df)
    for (cand in candidates){
      if (cand %in% nm) return(cand)
    }
    return(NA_character_)
  }
  
  sym_col <- map_col(df, col_symbol); if (is.na(sym_col)) stop("Could not find a symbol/Gene column.")
  tim_col <- map_col(df, col_time);   if (is.na(tim_col)) stop("Could not find a timepoint column.")
  lfc_col <- map_col(df, col_logfc);  if (is.na(lfc_col)) stop("Could not find a logFC column.")
  
  df <- df[, c(sym_col, tim_col, lfc_col)]
  names(df) <- c("symbol","timepoint","logFC")
  
  # Extract plain timepoint tokens (B0/D1/D5/F3/F6/F9/F12) from longer labels if needed
  df$timepoint <- as.character(df$timepoint)
  df$timepoint <- stringr::str_extract(df$timepoint, "(B0|D1|D5|F3|F6|F9|F12)")
  
  # Filter to requested proteins and timepoints
  df <- dplyr::filter(df, symbol %in% protein_names, timepoint %in% valid_times)
  
  if (nrow(df) == 0){
    stop("No matching rows for the requested proteins/timepoints.")
  }
  
  # Pivot wider to protein x time matrix
  mat <- tidyr::pivot_wider(df, id_cols = "symbol", names_from = "timepoint", values_from = "logFC")
  # Reorder columns by valid_times
  keep <- intersect(valid_times, names(mat))
  mat  <- mat[, c("symbol", keep), drop=FALSE]
  rownames(mat) <- mat$symbol
  mat$symbol <- NULL
  
  # Convert to numeric matrix
  M <- as.matrix(mat)
  storage.mode(M) <- "numeric"
  
  # Optionally z-score rows
  if (isTRUE(zscore_rows)){
    M <- t(scale(t(M)))
  }
  
  # Build a diverging palette centered at 0
  breaks <- 100
  # pheatmap centers colors at 0 if we pass appropriate colorRampPalette; we'll use default.
  
  # File path
  outfile <- file.path(output_dir, fname)
  
  # Plot and save
  pheatmap::pheatmap(
    M,
    filename = outfile,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(breaks),
    border_color = NA,
    main = "logFC (IF vs NIF)",
    na_col = "grey90",
    silent = TRUE
  )
  
  message("Saved heatmap: ", outfile)
  invisible(M)
}

#------------------------------ idk --------------------------------------

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

# -------------------------- Mark Idents -------------------------------------

outDir_Mark_idents <- "~/Roselab/Metabolite/results/Figures/MarkIdents/"

Mark_idents <- c("MW0004036|Isovanillic acid", "MW0115839|(2E)-3-(1H-indol-2-yl)prop-2-enoic acid",
                 "MW0017085|Cholesterol sulfate")

plot_metabolite_trends(metabolite_names = Mark_idents,
                       output_dir = outDir_Mark_idents,
                       metlData_log2 = metlData_log2,
                       metaData = metaData)

plot_metabolite_boxplots(metabolite_names = Mark_idents,
                         output_dir = outDir_Mark_idents,
                         metlData_log2 = metlData_log2,
                         metaData = metaData)


# -------------------------- Mark Idents -------------------------------------

outDir_Mark_idents <- "~/Roselab/Metabolite/results/Figures/MarkIdents/"

Mark_idents <- c("MW0004036|Isovanillic acid", "MW0115839|(2E)-3-(1H-indol-2-yl)prop-2-enoic acid",
                 "MW0017085|Cholesterol sulfate")

plot_metabolite_trends(metabolite_names = Mark_idents,
                       output_dir = outDir_Mark_idents,
                       metlData_log2 = metlData_log2,
                       metaData = metaData)

plot_metabolite_boxplots(metabolite_names = Mark_idents,
                         output_dir = outDir_Mark_idents,
                         metlData_log2 = metlData_log2,
                         metaData = metaData)
# -------------------------- Mark Idents2 -------------------------------------

outDir_Mark_idents2 <- "~/Roselab/Metabolite/results/Figures/MarkIdents/idents"

Mark_idents2 <- c("MW0147772|Coproporphyrin I",
                 "MW0010768|Cloprostenol",
                 "MW0126450|Riboflavin reduced",
                 "MW0105369|6-Sulfatoxymelatonin")

plot_metabolite_trends(metabolite_names = Mark_idents2,
                       output_dir = outDir_Mark_idents2,
                       metlData_log2 = metlData_log2,
                       metaData = metaData)

plot_metabolite_boxplots(metabolite_names = Mark_idents2,
                         output_dir = outDir_Mark_idents2,
                         metlData_log2 = metlData_log2,
                         metaData = metaData)



#-------------------------- Patient trengds ----------------------------

plot_metabolite_patient_trends(
  metabolite_names = "MW0126450|Riboflavin reduced",
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData,
  plot_title = "Patient-specific trajectories of Riboflavin reduced2"
)

plot_metabolite_patient_trends_overlap(
  metabolite_names = "MW0126450|Riboflavin reduced",
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData,
  plot_title = "Patient-specific trajectories of Riboflavin reduced"
)

"MW0155717|Prenylated FMNH2"

plot_metabolite_patient_trends_overlap(
  metabolite_names = "MW0155717|Prenylated FMNH2",
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData,
  plot_title = "Patient-specific trajectories of Prenylated FMNH2"
)

plot_metabolite_patient_trends_overlap(
  metabolite_names = "MW0155717|Prenylated FMNH2",
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData,
  plot_title = "Patient-specific trajectories of Prenylated FMNH2"
)


#------------------------- Rose Presentaion ---------------------------------
# c("B0", "D1", "D5", "F3", "F6")
Ferroptosis = c("MW0048971|Coenzyme Q10", "FDATN01488|L-Glutamic acid", 
                "MW0107690|L-cystine")

plot_metabolite_patient_trends_overlap(
  metabolite_names = Ferroptosis,
  valid_times = c("B0", "D1", "D5"),
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData
)


plot_metabolite_boxplots(
  metabolite_names = Ferroptosis,
  valid_times = c("B0", "D1", "D5"),
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData
)

long_study = c("MW0063752|Ursocholyltaurine", "MW0147772|Coproporphyrin I")
plot_metabolite_patient_trends_overlap(
  metabolite_names = long_study,
  valid_times = c("B0", "D1", "D5", "F3", "F6"),
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData
)

plot_metabolite_boxplots(
  metabolite_names = long_study,
  valid_times = c("B0", "D1", "D5", "F3", "F6"),
  output_dir = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  metlData_log2 = metlData_log2,
  metaData = metaData
)

#--------------------------- Heat maps ---------------------------------------
plot_heatmap <- function(
    protein_names = c("prot1","prot2"),      # genes/proteins to plot, order preserved on y-axis
    output_dir = ".",                         # output location (created if missing)
    valid_times = c("D1","D5","F3","F6","F9","F12"),  # timepoints of interest, order preserved on x-axis
    prot_results,                             # data.frame with columns: timepoint, Gene, logFC (others ignored)
    zscore_rows = FALSE,                      # z-score per gene across its timepoints
    fname = "heatmap_proteins.png"            # output filename (PNG)
) {
  # ---- dependencies (base + tidyverse) ----
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr', 'tidyr', and 'ggplot2' are required.")
  }
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---- basic checks ----
  needed <- c("timepoint","Gene","logFC")
  if (!all(needed %in% names(prot_results))) {
    stop("`prot_results` must contain columns: timepoint, Gene, logFC")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # keep only requested timepoints and proteins (preserve user order)
  tp_levels <- intersect(valid_times, unique(prot_results$timepoint))
  if (length(tp_levels) == 0) stop("None of the requested `valid_times` are present in `prot_results`.")
  prot_in_data <- intersect(protein_names, unique(prot_results$Gene))
  if (length(prot_in_data) == 0) stop("None of the requested `protein_names` are present in `prot_results`.")
  missing_prot <- setdiff(protein_names, prot_in_data)
  missing_tp   <- setdiff(valid_times, tp_levels)
  
  # summarize to one value per Gene × timepoint (mean if duplicates)
  long <- prot_results %>%
    filter(Gene %in% protein_names, timepoint %in% valid_times) %>%
    select(Gene, timepoint, logFC) %>%
    group_by(Gene, timepoint) %>%
    summarise(value = mean(logFC, na.rm = TRUE), .groups = "drop")
  
  # make complete grid so missing combos show as NA tiles
  long <- long %>%
    tidyr::complete(Gene = protein_names, timepoint = valid_times) %>%
    mutate(
      Gene = factor(Gene, levels = protein_names),
      timepoint = factor(timepoint, levels = valid_times)
    )
  
  # optional row-wise z-score
  if (zscore_rows) {
    long <- long %>%
      group_by(Gene) %>%
      mutate(value = if (all(is.na(value))) NA_real_ else
        as.numeric(scale(value))) %>%
      ungroup()
    fill_label <- "Row z-score"
  } else {
    fill_label <- "logFC"
  }
  
  # “pretty” diverging palette (blue→white→red), NA tiles light gray
  pretty_diverging <- grDevices::colorRampPalette(
    c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8",
      "#ffffbf",
      "#fee090","#fdae61","#f46d43","#d73027","#a50026")
  )
  
  p <- ggplot(long, aes(x = timepoint, y = Gene, fill = value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(colors = pretty_diverging(255), na.value = "grey90", name = fill_label) +
    labs(
      title = "logFC of Selected Proteins (IF vs NIF)",
      x = "Timepoint",
      y = "Protein"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 18, face = "bold"),
      axis.text.y   = element_text(size = 18, face = "bold"),
      axis.title.y  = element_text(size = 22, face = "bold"),
      axis.title.x  = element_text(size = 20, face = "bold"),
      plot.title    = element_text(hjust = 0.5, size = 25, face = "bold"),
      panel.grid    = element_blank(),
      legend.title  = element_text(size = 20, face = "bold"),   # bigger legend title
      legend.text   = element_text(size = 18),                  # bigger legend tick labels
      legend.key.height = unit(1.5, "cm"),                      # stretch legend bar
      legend.key.width  = unit(0.6, "cm")                       # make it wider
    )
  
  # dimensions scale with data size
  width_px  <- max(1200, 180 * length(valid_times))
  height_px <- max(900, 80 * length(protein_names))
  outfile <- file.path(output_dir, fname)
  ggplot2::ggsave(outfile, plot = p, width = width_px/96, height = height_px/96, dpi = 96, bg = "white")
  
  if (length(missing_prot) > 0 || length(missing_tp) > 0) {
    message(
      sprintf("Note: missing proteins: %s | missing timepoints: %s",
              ifelse(length(missing_prot)==0, "none", paste(missing_prot, collapse=", ")),
              ifelse(length(missing_tp)==0, "none", paste(missing_tp, collapse=", ")))
    )
  }
  
  invisible(list(plot_data = long, file = outfile))
}



# Protein Heatmap
my_prots <- c(
  "NPS-PLA2", # "PLA2G2A"
  "TFF1",
  "FCG2B",   #  FCGR2B
  "bFGF",     #  FGF2
  "PEDF",     # SERPINF1
  "LL-37",     #  CAMP
  "PDIA4",
  "WNT5A",
  #"AHSG",     # a2-HS-Glycoprotein
  "GSTK1",    #Glutathione S-transferase kappa 1
  "PSG3"       # Pregnancy zone protein
)
my_prots <- rev(my_prots)

plot_heatmap(
  protein_names = my_prots,
  output_dir    = "~/Roselab/Metabolite/results/Figures/MetabolitePlots/patient_line",
  valid_times   = c("D1","D5","F3","F6","F9","F12"),
  prot_results  = prot_results,
  fname = "Heatmap_prot.png"
)

#---------------------- Pulling specific data -------------------------------
a = metlData_log2["MW0155717|Prenylated FMNH2", ] 
b = data.frame(Sample_ID = colnames(metlData_log2), Value = as.numeric(a), 
               stringsAsFactors = FALSE) 
c = merge(b, metaData, by = "Sample_ID") 
c

write.csv(c, "~/Roselab/Metabolite/results/misc/Prenylated.csv")






