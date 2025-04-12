library(MetaboAnalystR)

cids <- get_cid("HMDB0013450", from = "hmdb")




run_enrichment_analysis <- function(results, output_dir, which_results, direction = NULL) {
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  # Load MetaboAnalystR if not already loaded
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    stop("Package 'MetaboAnalystR' is required but not installed.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Filter your results based on the logFC direction parameter
  if (is.null(direction)) {
    hmdb_list <- results %>%
      filter(P.Value <= 0.05 & HMDB != "-") %>%
      pull(HMDB)
    plt_name <- "ALL"
  } else if (toupper(direction) == "UP") {
    hmdb_list <- results %>%
      filter(P.Value <= 0.05 & log2FC > 0 & HMDB != "-") %>%
      pull(HMDB)
    plt_name <- "UP"
  } else if (toupper(direction) == "DOWN") {
    hmdb_list <- results %>%
      filter(P.Value <= 0.05 & log2FC < 0 & HMDB != "-") %>%
      pull(HMDB)
    plt_name <- "DOWN"
  } else {
    stop("Invalid 'direction' parameter provided. Use NULL, 'UP', or 'DOWN'.")
  }
  
  if (length(hmdb_list) == 0) {
    warning("No metabolites passed the filtering criteria. Returning NULL.")
    return(NULL)
  }
  
  # Step 2: Initialize MetaboAnalyst object for enrichment analysis
  mSet <- MetaboAnalystR::InitDataObjects("conc", "msetora", FALSE)
  
  # Step 3: Map the filtered compound list
  mSet <- MetaboAnalystR::Setup.MapData(mSet, hmdb_list)
  
  # Step 4: Cross-reference the compounds using HMDB IDs
  mSet <- MetaboAnalystR::CrossReferencing(mSet, "hmdb")
  
  # Step 5: Generate a mapping table to see which metabolites matched to which pathway
  mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
  
  # Optional: if you want to check/correct mappings manually, you could add additional calls here
  
  # Step 6: Disable metabolome filter (turn off background trimming)
  mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
  
  print("error here?")
  print(mSet)
  # Step 7: Set the metabolite set library for pathway analysis (e.g., SMPDB pathways)
  mSet <- MetaboAnalystR::SetCurrentMsetLib(mSet, "smpdb_pathway", 0)
  print("before ehre?")
  
  # Defensive checks before CalculateHyperScore()
  cat("Checking compound-pathway mapping before enrichment...\n")
  
  # Check mapping table
  if (is.null(mSet$dataSet$map.table) || nrow(mSet$dataSet$map.table) == 0) {
    stop("Mapping table is empty. No compounds were successfully mapped.")
  }
  
  
  # Step 8: Run the enrichment analysis (hypergeometric test)
  mSet <- MetaboAnalystR::CalculateHyperScore(mSet)
  
  print("calculated hyperscore")
  
  # Check if enrichment results exist and are properly formatted
  enrichment_results <- mSet$analSet$ora.mat
  if (is.null(enrichment_results) || !is.matrix(enrichment_results)) {
    warning("No enrichment results found or result is not a matrix. Possibly no hits.")
    return(NULL)
  }
  
  # Extract the enrichment analysis results table
  enrichment_results <- mSet$analSet$ora.mat
  
  # Prepare data frame for plotting (include pathway names)
  ora_df <- as.data.frame(enrichment_results)
  ora_df$Pathway <- rownames(ora_df)
  
  ## ---------------------- Bar Plot ----------------------
  # Filter out pathways with at least one hit and order by Raw p-value
  bar_df <- ora_df %>%
    filter(hits > 0) %>%
    arrange(`Raw p`)
  top_bar_df <- head(bar_df, 10)  # Top 10 pathways
  
  # Create bar plot showing -log10(p-value)
  bar_plot <- ggplot(top_bar_df, aes(x = reorder(Pathway, -`Raw p`), y = -log10(`Raw p`))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Pathway", y = "-log10(p-value)", title = "Top Enriched Pathways") +
    theme_bw()
  
  # Save the bar plot to the specified output directory
  bar_plot_file <- file.path(output_dir, paste0(which_results, "_ORA_barplot", plt_name, ".png"))
  ggsave(filename = bar_plot_file, plot = bar_plot, dpi = 300)
  
  ## ---------------------- Dot Plot ----------------------
  # Add additional columns for the dot plot
  ora_df <- ora_df %>%
    mutate(log_p = -log10(`Raw p`),
           hit_ratio = hits / total) %>%
    filter(hits > 0) %>%
    arrange(`Raw p`) %>%
    slice(1:15)  # top 15 pathways
  
  # Create dot plot of enrichment results
  dot_plot <- ggplot(ora_df, aes(x = log_p, y = reorder(Pathway, log_p))) +
    geom_point(aes(size = hit_ratio), color = "steelblue") +
    scale_size_continuous(name = "Hit Ratio (hits/total)") +
    labs(
      x = expression(-log[10](p)),
      y = "Pathway",
      title = "Enriched Pathways (ORA)"
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10))
  
  # Save the dot plot
  dot_plot_file <- file.path(output_dir, paste0(which_results, "_ORA_dotplot", plt_name, ".png"))
  ggsave(filename = dot_plot_file, plot = dot_plot, dpi = 300)
  
  # Return a list containing the analysis object, enrichment results,
  # and the file paths to the saved plots. You can further inspect the mSet object
  # to see which metabolites matched to which pathways.
  return(list(
    mSet = mSet,
    enrichment_results = enrichment_results,
    bar_plot_file = bar_plot_file,
    dot_plot_file = dot_plot_file
  ))
}


# # test_extract <- read.csv("~/Roselab/Metabolite/results/MarkResults/between_group_pvalues_F3.csv")
# my_results <- read.csv("~/Roselab/Metabolite/results/MarkResults/between_group_pvalues_F3.csv")
# # my_results <- read.csv("~/Roselab/Metabolite/results/difference/limma/intraArm/limma_intraArm_F3_IF_vs_NIF.csv")
# 
# output_fig <- "~/Roselab/Metabolite/results/Figures/ORA/"
# 
# analysis_output <- run_enrichment_analysis(results = my_results,
#                                            output_dir = output_fig,
#                                            which_results = "F3",
#                                            direction = "UP")
# # Check the enrichment table:
# head(analysis_output$enrichment_results)
# 
# 
# analysis_output$mSet$analSet$ora.hits$`Alpha Linolenic Acid and Linoleic Acid Metabolism`
# 
# mapped_table <- analysis_output$mSet$dataSet$map.table
# mapped_table[mapped_table[, "Match"] == "Cis-8,11,14,17-Eicosatetraenoic acid", ]



#------------------- Common metabolites --------------------------



# test_extract <- read.csv("~/Roselab/Metabolite/results/MarkResults/between_group_pvalues_F3.csv")
my_results_common <- read.csv("~/Roselab/Metabolite/results/difference/common_D5_F3/D5_common_metabolite_results.csv")
# my_results <- read.csv("~/Roselab/Metabolite/results/difference/limma/intraArm/limma_intraArm_F3_IF_vs_NIF.csv")

output_fig_common <- "~/Roselab/Metabolite/results/Figures/MetabolitePlots/common_metabolites/ORA/"

analysis_output <- run_enrichment_analysis(results = my_results_common,
                                           output_dir = output_fig_common,
                                           which_results = "Common_metabolites"
                                           )
# Check the enrichment table:
head(analysis_output$enrichment_results)






















#------------------ Test_data ------------------------------------

# Step 1: Prepare your vector of HMDB IDs (already done)
hmdb_list <- test_extract %>%
  filter(p_value <= 0.05 & HMDB != "-") %>%
  pull(HMDB)
length(hmdb_list)

hmdb_list <- test_extract %>%
  filter(p_value <= 0.05 & log2FC > 0 & HMDB != "-",) %>%
  pull(HMDB)
length(hmdb_list)

hmdb_list <- test_extract %>%
  filter(p_value <= 0.05 & log2FC < 0 & HMDB != "-",) %>%
  pull(HMDB)
length(hmdb_list)



# Step 2: Initialize MetaboAnalyst object for enrichment
mSet <- InitDataObjects("conc", "msetora", FALSE)

# Step 3: Map the compound list
mSet <- Setup.MapData(mSet, hmdb_list)

# Step 4: Cross-reference using HMDB IDs
mSet <- CrossReferencing(mSet, "hmdb")

# Step 5: Create mapping table
mSet <- CreateMappingResultTable(mSet)

# Optional: Check or correct mappings manually
# e.g. mSet <- PerformDetailMatch(mSet, "HMDB0001234")
# mSet <- GetCandidateList(mSet)
# mSet <- SetCandidate(mSet, "HMDB0001234", "Correct_Name")

# Step 6: Disable metabolome filter (optional, disables background trimming)
mSet <- SetMetabolomeFilter(mSet, FALSE)

# Step 7: Choose metabolite set library
# Options include: "smpdb_pathway", "kegg_pathway", "smpdb_disease", etc.
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 0)

# Step 8: Run enrichment
mSet <- CalculateHyperScore(mSet)

# Step 9: View results
head(mSet$analSet$ora.mat)

# Step 10: Optional - Plot ORA
mSet <- PlotORA(mSet, "ora_plot", "bar", "png", 72)


library(ggplot2)
ora_df <- as.data.frame(mSet$analSet$ora.mat)
ora_df$Pathway <- rownames(ora_df)

# Filter for hits and sort
ora_df <- ora_df[ora_df$hits > 0, ]
ora_df <- ora_df[order(ora_df$`Raw p`), ]

# Plot top 10
ggplot(head(ora_df, 10), aes(x = reorder(Pathway, -`Raw p`), y = -log10(`Raw p`))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Pathway", y = "-log10(p-value)", title = "Top Enriched Pathways") +
  theme_minimal()



# ----------------------------- Dot plot -----------------------------------

# Extract and format enrichment results
ora_df <- as.data.frame(mSet$analSet$ora.mat)
ora_df$Pathway <- rownames(ora_df)

# Create new columns for plotting
ora_df <- ora_df %>%
  mutate(log_p = -log10(`Raw p`),
         hit_ratio = hits / total)

# Optional: filter top N (e.g., top 15 by p-value)
ora_plot_df <- ora_df %>%
  filter(hits > 0) %>%
  arrange(`Raw p`) %>%
  slice(1:15)

# Plot
ggplot(ora_plot_df, aes(x = log_p, y = reorder(Pathway, log_p))) +
  geom_point(aes(size = hit_ratio), color = "steelblue") +
  scale_size_continuous(name = "Hit Ratio (hits/total)") +
  labs(
    x = expression(-log[10](p)),
    y = "Pathway",
    title = "Enriched Pathways (ORA)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

