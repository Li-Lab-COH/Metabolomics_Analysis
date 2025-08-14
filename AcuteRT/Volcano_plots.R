library(ggplot2)

make_volcano_plots <- function(
    in_dir  = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/LimmaResults/AllResults/",
    out_dir = "~/1Work/Roselab/Metabolomics/results/acute_vs_chronic_RT/figs/volcano",
    p_col_candidates   = c("P.Value", "adj.P.Val"),
    lfc_col = "logFC",
    p_cut   = 0.05,
    lfc_cut = 1,     # log2 fold-change threshold
    width   = 8,
    height  = 6,
    dpi     = 300
){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("No CSV files found in: ", in_dir)
    return(invisible(NULL))
  }
  
  for (f in files) {
    df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(df)) { warning("Could not read: ", f); next }
    
    # Find columns
    if (!lfc_col %in% colnames(df)) { warning("Missing logFC in ", f); next }
    p_col <- p_col_candidates[p_col_candidates %in% colnames(df)][1]
    if (is.na(p_col)) { warning("No P.Value/adj.P.Val in ", f); next }
    
    # Coerce to numeric just in case
    df[[lfc_col]] <- as.numeric(df[[lfc_col]])
    df[[p_col]]   <- as.numeric(df[[p_col]])
    df$negLog10P  <- -log10(df[[p_col]])
    
    # Significance flag
    df$signif <- with(df, ifelse(abs(df[[lfc_col]]) >= lfc_cut & df[[p_col]] <= p_cut, "sig", "ns"))
    
    # Build plot
    p <- ggplot(df, aes(x = .data[[lfc_col]], y = negLog10P)) +
      geom_point(aes(color = signif), alpha = 0.7, size = 1.8) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
      scale_color_manual(values = c("ns" = "grey70", "sig" = "firebrick")) +
      labs(
        title = tools::file_path_sans_ext(basename(f)),
        x = "log2 fold change",
        y = "-log10(p)"
      ) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
    
    out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), "_volcano.png"))
    ggsave(out_file, plot = p, width = width, height = height, dpi = dpi)
    message("Saved: ", out_file)
  }
}

# Run it
make_volcano_plots()
