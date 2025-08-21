# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(writexl)
})



# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12

# =============================================================================
# Loading data

#Loading Metabolite data
# metlData <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_data_74.csv", row.names = "ID")
metlData <- read.csv(
  "~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_data_74.csv",
  row.names = "ID",
  check.names = FALSE
)
metl_metaData <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/meta_data_74.csv")
metabolite_mol_meta_data <- read.csv("~/1Work/Roselab/Metabolomics/data/metabolomics/metabolite_metadata.csv")
rownames(metlData) <- sub("\\|.*", "", rownames(metlData))


# Adding patient data
metl_metaData <- metl_metaData %>%
  mutate(Patient_ID = str_extract(Sample_ID, "^IF\\d+")) %>%
  filter(!is.na(Patient_ID) & Patient_ID !="")

# Metabolites of interest
# Cer(d18:0/12:0) - MW0055322
# LPC(22:5/0:0) - MEDP1322
# PE-NMe(18:1(9Z)/15:0) - MW0059344
# PE-NMe2(16:0/18:1(11Z)) - MW0060249


# =============================================================================
# Pseudocode

# Dataprep function
# 1. For each patient that has data in D1 and B0 use those
  # So if there is a patient that has metadata for both D1 and B0 in the metadata, 
  # grab those sample names and grab data from metalData for those sample names (columns).
  # These are the sample names that we will be using to calculate the
  # ratio of change. For example, in the metadata that I'm sharing IF05 has B0 as the sample ID: IF05_B1
  # and IF05_B2_P2 as the sampleID for D1. Eventually I want to calculate the ratio of change as IF05_B1/IF05_B2_P2.
  # The index of this dataframe should have the patient ID name, or should it be a column where we define which row
  # belongs to which patient? Considering that I want to isolate the datapoints for statistical comparison, should we 
  # make a column of patient ID and what cohort this patient belongs to, to make things easier in the future?
# 2. Isolate the B0 and D1 data for the metabolites of interest. 
  # This should be customizable in the function, if I wanted to I should be able to 
  # grab all of the timepoints
# 3. For each metabolite create a dataframe of the form:
  # * Each row represents one unique patient ID
  # * Column 1, will contain data for B0
  # * Column 2, will contain data for D1
  # * Column n, will contain data for [timpoits of interest]
  # * Final columns will be a calculation of the ratio of change (D1/B0). defined in the function. 
  #   The ratio of change of interest should be able to be defined in the function. For example, 
  #   c(D1, B0) would mean D1/B0 and the column would be named D1_B0
  #   c(D5, B0) would mean D5/B0 and the column would be named D5_B0
# 4. Output a dataframe for each metabolite of interest and have the dataframe
  # be named after the index given. I am flexible on the structure of data that is outputted.
  # The only requirement that it should be accessbile for plotting the ratio of change and
  # be able to isolate the patients that are either NIF or IF

# Statistics
# 1. Using one of the datasets as input from the previous function
  # extract the column of interest (e.g. D1_B0 - user inputed) and separate the data into two groups:
  # NIF and IF. We can use the metadata for this or column (see text in point 1 of dataprep function)
# 2. Use those two groups to do a pair wise t-test
# 3. output results


# Plotting
# 1. I will want to plot the ratio of change D1_B0 and D5_B0, or what ever column that I want plot
  # This has be an input of the function
# 2. Where for each column we group the data for each cohort, NIF or IF, plot the points, connect each 
  # column by the mean (D1_B0 to D5_B0) and the plot the standard error mean bars


# =============================================================================
# Functions
Isolate_data_function <- function(
    data,
    metaData,
    index
    ){
  data_log <- log2(data + 1)
  testing <- data_log[index, c("IF05_B1", "IF05_B2_P2")]
  return(testing)
}


Isolate_data_function(
  metlData,
  metl_metaData,
  index = "MW0055322"
)



# =============================================================================
# Scratch coding
stopifnot(
  all(metl_metaData %>% group_by(Patient_ID) %>% summarise(n = n_distinct(Arm)) %>% pull(n) == 1)
)

metlData_log = log2(metlData + 1)

# Inputs you control
times_of_interest  <- c("B0","D1")
index_of_interest  <- "MW0055322"

# 1) Keep only patients that have all requested times
map <- metl_metaData %>%
  filter(Time %in% times_of_interest) %>%
  group_by(Patient_ID) %>%
  filter(n_distinct(Time) == length(times_of_interest)) %>%  # require all times present
  ungroup() %>%
  select(Patient_ID, Arm, Time, Sample_ID) # ADD COHORT INFO

# Pull the values for the index from data_log, one per Sample_ID
stopifnot(index_of_interest %in% rownames(metlData_log))
map$Value <- as.numeric(metlData_log[index_of_interest, map$Sample_ID, drop = TRUE])


out <- map %>%
  group_by(Patient_ID, Time) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    Arm   = dplyr::first(Arm),   # Arm constant within a patient
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols    = c(Patient_ID, Arm),   # Arm stays as a column
    names_from = Time,
    values_from = Value
  ) %>%
  arrange(Patient_ID)




# =============================================================================
# Testing Functions



#-----------------------------
# Helpers
#-----------------------------
.sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

#-----------------------------
# 1) Build per-patient matrix for selected times
#-----------------------------
build_patient_time_matrix <- function(expr_data,
                                      metaData,
                                      times_of_interest = c("B0","D1"),
                                      sample_col = "Sample_ID",
                                      time_col   = "Time",
                                      arm_col    = "Arm",
                                      patient_col = "Patient_ID",
                                      transform_y = c("none","log2p1")) {
  transform_y <- match.arg(transform_y)
  
  # numeric sample columns only
  num_cols <- vapply(expr_data, is.numeric, logical(1))
  X <- as.matrix(expr_data[, num_cols, drop = FALSE])
  
  # (optional) transform
  if (transform_y == "log2p1") X <- log2(X + 1)
  
  # keep only rows of metadata whose samples exist in X
  meta_ok <- metaData %>% filter(.data[[sample_col]] %in% colnames(X))
  # filter to times of interest
  meta_ok <- meta_ok %>% filter(.data[[time_col]] %in% times_of_interest)
  
  # require patients to have ALL requested times
  keep_patients <- meta_ok %>%
    group_by(.data[[patient_col]]) %>%
    summarise(n_times = n_distinct(.data[[time_col]]), .groups="drop") %>%
    filter(n_times == length(times_of_interest)) %>%
    pull(.data[[patient_col]])
  
  meta_keep <- meta_ok %>% filter(.data[[patient_col]] %in% keep_patients)
  
  # Map (Patient, Arm, Time, Sample_ID)
  map <- meta_keep %>%
    select(Patient_ID = all_of(patient_col),
           Arm        = all_of(arm_col),
           Time       = all_of(time_col),
           Sample_ID  = all_of(sample_col))
  
  # Pull values from X for *each* Sample_ID, for each requested metabolite later.
  # Here we just return the map and X; extraction happens downstream to avoid
  # copying matrices repeatedly.
  list(X = X, map = map, times = times_of_interest)
}

#-----------------------------
# 2) Extract one metabolite’s values, average replicates, reshape wide,
#    and add ratio columns
#-----------------------------
prep_one_metabolite <- function(index_id,
                                build_obj,  # output of build_patient_time_matrix
                                ratios = list(c("D1","B0")),  # list of 2-length character vectors
                                aggregation = mean) {
  
  X   <- build_obj$X
  map <- build_obj$map
  times_of_interest <- build_obj$times
  
  stopifnot(index_id %in% rownames(X))
  
  # Values for each Sample_ID
  vals <- as.numeric(X[index_id, map$Sample_ID, drop = TRUE])
  map$Value <- vals
  
  # Average replicates within (Patient, Time); carry Arm
  long <- map %>%
    group_by(Patient_ID, Arm, Time) %>%
    summarise(Value = aggregation(Value, na.rm = TRUE), .groups = "drop")
  
  # Wide: Patient x Time
  wide <- long %>%
    pivot_wider(id_cols = c(Patient_ID, Arm),
                names_from = Time,
                values_from = Value)
  
  # Add ratio columns (e.g., D1/B0) with safe division
  for (pair in ratios) {
    stopifnot(length(pair) == 2)
    num <- pair[1]; den <- pair[2]
    new_name <- paste0(num, "_", den)
    if (!all(c(num, den) %in% colnames(wide))) {
      warning("Skipping ratio ", new_name, " because one of [", num, ", ", den, "] is missing.")
      next
    }
    wide[[new_name]] <- wide[[num]] / wide[[den]]
  }
  
  # Keep consistent ordering
  wide <- wide %>% arrange(Patient_ID)
  wide
}

#-----------------------------
# Test
#-----------------------------
compare_ratio_groups <- function(df_ratio,
                                 ratio_col,
                                 arm_col = "Arm",
                                 arms = c("NIF","IF"),
                                 method = c("t","ranksum","both"),
                                 paired = FALSE,
                                 var_equal = FALSE,        # for Student's t; FALSE = Welch
                                 alternative = "two.sided",
                                 conf_level = 0.95) {
  method <- match.arg(method)
  
  d <- df_ratio %>% dplyr::filter(.data[[arm_col]] %in% arms)
  g1 <- d %>% dplyr::filter(.data[[arm_col]] == arms[1]) %>% dplyr::pull(dplyr::all_of(ratio_col))
  g2 <- d %>% dplyr::filter(.data[[arm_col]] == arms[2]) %>% dplyr::pull(dplyr::all_of(ratio_col))
  
  out <- list()
  
  if (method %in% c("t","both")) {
    tt <- t.test(g1, g2, paired = paired, var.equal = var_equal,
                 alternative = alternative, conf.level = conf_level)
    out$t <- tibble::tibble(
      method = if (paired) "t_paired" else if (var_equal) "t_student" : "t_welch",
      ratio_col = ratio_col,
      group1 = arms[1], group2 = arms[2],
      n1 = sum(!is.na(g1)), n2 = sum(!is.na(g2)),
      mean1 = mean(g1, na.rm = TRUE), mean2 = mean(g2, na.rm = TRUE),
      diff  = mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE),
      statistic = unname(tt$statistic),
      parameter = unname(tt$parameter),
      p.value   = unname(tt$p.value),
      estimate  = if (!is.null(tt$estimate)) unname(diff(tt$estimate)) else NA_real_,
      conf.low  = if (!is.null(tt$conf.int)) unname(tt$conf.int[1]) else NA_real_,
      conf.high = if (!is.null(tt$conf.int)) unname(tt$conf.int[2]) else NA_real_,
      alternative = alternative
    )
  }
  
  if (method %in% c("ranksum","both")) {
    wx <- wilcox.test(g1, g2,
                      paired = paired,                # FALSE = rank-sum; TRUE = signed-rank
                      alternative = alternative,
                      conf.int = TRUE, conf.level = conf_level,
                      exact = FALSE, correct = TRUE)  # robust to ties/large n
    out$w <- tibble::tibble(
      method = if (paired) "wilcoxon_signed_rank" else "wilcoxon_rank_sum",
      ratio_col = ratio_col,
      group1 = arms[1], group2 = arms[2],
      n1 = sum(!is.na(g1)), n2 = sum(!is.na(g2)),
      median1 = stats::median(g1, na.rm = TRUE),
      median2 = stats::median(g2, na.rm = TRUE),
      hl_estimate = if (!is.null(wx$estimate)) unname(wx$estimate) else NA_real_,  # Hodges–Lehmann
      statistic = unname(wx$statistic),
      p.value   = unname(wx$p.value),
      conf.low  = if (!is.null(wx$conf.int)) unname(wx$conf.int[1]) else NA_real_,
      conf.high = if (!is.null(wx$conf.int)) unname(wx$conf.int[2]) else NA_real_,
      alternative = alternative
    )
  }
  
  dplyr::bind_rows(out)
}


#-----------------------------
# 4) Plot one or more ratio columns, grouped by Arm,
#    points + mean line + SEM error bars, fixed dimensions.
#-----------------------------
plot_ratio_columns <- function(df_ratio,
                               ratio_cols = c("D1_B0","D5_B0"),
                               arm_col = "Arm",
                               title_main = "",
                               subtitle = "",
                               index_id,
                               output_dir = "ratio_plots",
                               fig_width_mm = 70,
                               fig_height_mm = 50,
                               base_pt = 7,
                               jitter_width = 0.12,   # horizontal jitter
                               jitter_height = 0,     # vertical jitter
                               dodge_width = 0.45) {  # distance between cohorts
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (length(ratio_cols) == 0) stop("ratio_cols is empty.")
  
  long <- df_ratio %>%
    dplyr::select(Patient_ID, dplyr::all_of(arm_col), dplyr::all_of(ratio_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(ratio_cols), names_to = "Ratio", values_to = "Value") %>%
    dplyr::mutate(Ratio = factor(Ratio, levels = ratio_cols))
  
  sumdf <- long %>%
    dplyr::group_by(Ratio, .data[[arm_col]]) %>%
    dplyr::summarise(meanValue = mean(Value, na.rm = TRUE),
                     semValue  = stats::sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
                     n = sum(!is.na(Value)),
                     .groups = "drop")
  
  # sizing knobs (mm for lines; point size is in mm-ish units)
  pt_mean_mm  <- 1.3
  pt_jitter_mm <- 1.0
  line_lw_mm  <- 0.3
  eb_cap_mm   <- 0.5
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = Ratio, y = Value, color = .data[[arm_col]])) +
    # jittered individual points per cohort (so overlapping points separate)
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = jitter_width,
                                               jitter.height = jitter_height,
                                               dodge.width = dodge_width),
      alpha = 0.8, size = pt_jitter_mm
    ) +
    # mean point per cohort
    ggplot2::geom_point(
      data = sumdf,
      ggplot2::aes(x = Ratio, y = meanValue, color = .data[[arm_col]], group = .data[[arm_col]]),
      position = ggplot2::position_dodge(width = dodge_width),
      shape = 18, size = pt_mean_mm, inherit.aes = FALSE
    ) +
    # SEM bars, colored by cohort
    ggplot2::geom_errorbar(
      data = sumdf,
      ggplot2::aes(x = Ratio, ymin = meanValue - semValue, ymax = meanValue + semValue,
                   color = .data[[arm_col]], group = .data[[arm_col]]),
      width = eb_cap_mm, linewidth = line_lw_mm,
      position = ggplot2::position_dodge(width = dodge_width),
      inherit.aes = FALSE
    ) +
    # (mean connecting line REMOVED)
    ggplot2::labs(title = title_main, subtitle = subtitle,
                  x = "Ratio", y = "Fold-change (unitless)", color = arm_col) +
    ggplot2::theme_bw(base_size = base_pt)
  
  file_base <- file.path(output_dir, paste0(index_id, "_", paste(ratio_cols, collapse = "_")))
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot = p, width = fig_width_mm, height = fig_height_mm, units = "mm")
  ggplot2::ggsave(paste0(file_base, ".png"), plot = p, width = fig_width_mm, height = fig_height_mm, units = "mm", dpi = 300)
  
  invisible(p)
}

#-----------------------------
# 5) Convenience wrapper for multiple metabolites
#-----------------------------
process_and_plot_metabolites <- function(index_ids,
                                         expr_data,
                                         metaData,
                                         molMetadata = NULL,
                                         times = c("B0","D1","D5"),
                                         ratios = list(c("D1","B0"), c("D5","B0")),
                                         transform_y = "log2p1",
                                         output_dir = "ratio_plots",
                                         stat_method = c("t","ranksum","both"),
                                         paired = FALSE,
                                         var_equal = FALSE,
                                         alternative = "two.sided") {
  stat_method <- match.arg(stat_method)
  
  builder <- build_patient_time_matrix(expr_data = expr_data,
                                       metaData  = metaData,
                                       times_of_interest = times,
                                       transform_y = transform_y)
  
  # Optional: Index -> pretty compound name
  title_lookup <- NULL
  if (!is.null(molMetadata) && all(c("Index","Compounds") %in% colnames(molMetadata))) {
    title_lookup <- setNames(molMetadata$Compounds, molMetadata$Index)
  }
  
  results <- list()
  
  for (idx in index_ids) {
    df_idx <- prep_one_metabolite(index_id = idx, build_obj = builder, ratios = ratios)
    
    # figure which ratio columns exist
    ratio_cols_present <- setdiff(colnames(df_idx), c("Patient_ID","Arm", times))
    
    tt_out <- dplyr::bind_rows(lapply(
      ratio_cols_present,
      function(rc) compare_ratio_groups(df_idx, ratio_col = rc,
                                        method = stat_method,
                                        paired = paired,
                                        var_equal = var_equal,
                                        alternative = alternative)
    ))
    
    # titles: pretty name in title, index as subtitle
    pretty_title <- if (!is.null(title_lookup) && idx %in% names(title_lookup)) title_lookup[[idx]] else idx
    
    # ✅ pass idx for safe filenames
    plot_ratio_columns(df_ratio = df_idx,
                       ratio_cols = ratio_cols_present,
                       title_main = pretty_title,
                       subtitle  = idx,
                       index_id  = idx,
                       output_dir = output_dir)
    
    results[[idx]] <- list(data = df_idx, ttest = tt_out)
  }
  
  results
}
save_ratio_results_excel <- function(results,
                                     index_id,
                                     output_dir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/ratio_tests/",
                                     add_test_col = TRUE) {
  # require writexl (install if needed)
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required. Install with install.packages('writexl').")
  }
  if (!index_id %in% names(results)) stop("index_id not found in results.")
  
  df <- results[[index_id]]$ttest
  if (is.null(df)) stop("No $ttest table found for index: ", index_id)
  
  # Add your test label (Wilcoxon vs Welch) if requested
  if (add_test_col && !"test" %in% names(df)) {
    df <- dplyr::mutate(df, test = ifelse(is.na(.data$parameter), "wilcoxon_rank_sum", "t_welch"))
  }
  
  # Ensure directory exists
  out_dir <- path.expand(output_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Save as Excel named after the index
  out_file <- file.path(out_dir, paste0(index_id, "_ratio_tests.xlsx"))
  writexl::write_xlsx(df, out_file)
  
  message("Saved: ", out_file)
  invisible(out_file)
}



#========================================
# EXAMPLE USAGE FOR YOUR THREE METABOLITES
#========================================

# Define metabolites of interest (Index IDs)
index_ids <- c("MW0055322", "MW0054553", "MW0059439")

# Choose times and ratios
times_to_use <- c("B0","D1","D5")
ratio_list   <- list(c("D1","B0"), c("D5","B0"))

out_figs_metl_ratio <- "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/figs/metl/ratio"

# Run (uses log2p1 by default; set to "none" to use raw scale)
# Assumes metlData, metl_metaData already in memory;
# optional molMetadata with columns: Index, Compounds for titles
res <- process_and_plot_metabolites(
  index_ids   = c("MW0055322","MEDP1322","MW0059344", "MW0060249"),
  expr_data   = metlData,
  metaData    = metl_metaData,
  molMetadata = metabolite_mol_meta_data,
  times       = c("B0","D1","D5"),
  ratios      = list(c("D1","B0"), c("D5","B0")),
  transform_y = "log2p1",
  output_dir  = out_figs_metl_ratio,
  stat_method = "both" #"ranksum"           # <— non-parametric rank-sum
)


# Access outputs:
# res[["MW0055322"]]$data    # per-patient wide table with B0, D1, D5, D1_B0, D5_B0
# res[["MW0055322"]]$ttest   # t-test summary for D1_B0 and D5_B0
"~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/ratio_tests/"

save_ratio_results_excel(res, "MW0055322")
save_ratio_results_excel(res, "MEDP1322")
save_ratio_results_excel(res, "MW0059344")
save_ratio_results_excel(res, "MW0060249")
