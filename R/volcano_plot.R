volcano_plot <- function(X,
                         showFDRLine = TRUE,
                         alpha = 0.05) {

  # Calculate FDR-adjusted p-values for determining the adjusted p-value line
  X_adjusted <- dplyr::mutate(X, p_adj = stats::p.adjust(pvalue, method = "fdr"))

  # Determine the raw p-value that corresponds to an adjusted p-value of alpha
  # This requires finding the point where p_adj is closest to alpha.
  # We sort by p_adj and then find the first pvalue where p_adj >= alpha
  fdr_pvalue_threshold <- NA
  if (showFDRLine == TRUE) {
    # Sort the data by adjusted p-value to find the threshold
    sorted_X_adjusted <- dplyr::arrange(X_adjusted, p_adj)
    # Find the row where p_adj first crosses or equals alpha
    # If no p_adj is less than or equal to alpha, this will be NA
    fdr_threshold_row <- dplyr::filter(sorted_X_adjusted, p_adj <= alpha)
    if (nrow(fdr_threshold_row) > 0) {
      # Take the maximum raw p-value among those that are significant after FDR adjustment
      fdr_pvalue_threshold <- max(fdr_threshold_row$pvalue)
    }
  }

  # Use FDR-adjusted or raw p-values as the significance measure for coloring points
  if (showFDRLine == TRUE) {
    X <- dplyr::mutate(X, p_sig = stats::p.adjust(pvalue, method = "fdr"))
  } else {
    X <- dplyr::mutate(X, p_sig = pvalue)
  }

  # Assign significance labels based on the threshold
  X <- dplyr::mutate(
    X,
    Significant = dplyr::case_when(
      p_sig > alpha ~ "NS",
      p_sig <= alpha & log2FoldChange >= 0 ~ "Up",
      p_sig <= alpha & log2FoldChange < 0 ~ "Down"
    )
  )
  X$Significant <- factor(X$Significant, levels = c("Down", "NS", "Up"))

  # Create labels with counts for each significance group
  # These are now calculated *after* X$Significant has been assigned
  ns_count <- sum(X$Significant == "NS")
  up_count <- sum(X$Significant == "Up")
  down_count <- sum(X$Significant == "Down")

  # Define the labels for the legend explicitly
  legend_labels <- c(
    "Down" = paste0("Down (", down_count, ")"),
    "NS" = paste0("NS (", ns_count, ")"),
    "Up" = paste0("Up (", up_count, ")")
  )

  # Define the colors for the legend explicitly
  legend_colors <- c("Down" = "blue", "NS" = "black", "Up" = "red")


  # Select the top 10 points with the highest absolute log2 fold change among significant genes

  top10 <- dplyr::pull(dplyr::slice_head(dplyr::arrange(
    dplyr::filter(X, p_sig <= alpha), dplyr::desc(abs(log2FoldChange))
  ), n = 10), Symbol)

  # Label top 10 points for the plot

  X <- dplyr::mutate(X,
                     top10label = dplyr::if_else(Symbol %in% top10, Symbol, NA_character_))

  # Build the volcano plot

  p <- ggplot2::ggplot(X,
                       ggplot2::aes(
                         x = log2FoldChange,
                         y = -log10(pvalue),
                         col = Significant,
                         label = top10label
                       )) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(labels = legend_labels,
                                 values = legend_colors) +
    ggrepel::geom_text_repel(ggplot2::aes(fontface = "bold"),
                             size = 9 / ggplot2::.pt,
                             na.rm = TRUE) +
    ggplot2::geom_hline(
      linetype = "dotted",
      yintercept = -log10(alpha), # Line for raw p-value 0.05
      col = "black",
      key_glyph = "blank" # Hide from legend
    ) +
    ggprism::theme_prism(base_size = 9) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.spacing.y = grid::unit(-1, "mm"),
      legend.spacing.x = grid::unit(-1, "mm"),
      legend.box.spacing = grid::unit(-1, "mm"),
      legend.key.spacing.x = grid::unit(-1, "mm"),
      legend.key.spacing.y = grid::unit(-1, "mm"),
      legend.box.margin = ggplot2::margin(),
      legend.box = "vertical",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 1.5, 0, 0), "mm"),
      panel.spacing = grid::unit(c(0, 0, 0, 0), "mm"),
      legend.margin = ggplot2::margin(0),
      legend.text = ggplot2::element_text(
        face = "bold",
        size = 9,
        margin = ggplot2::margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0
        )
      )
    )

  # If using FDR, add an additional horizontal line at the p-value corresponding to the FDR threshold
  if (showFDRLine == TRUE && !is.na(fdr_pvalue_threshold)) {
    p <- p + ggplot2::geom_hline(
      linetype = "dashed", # Use a different linetype to distinguish
      yintercept = -log10(fdr_pvalue_threshold),
      col = "black" # Use a different color for the FDR line
    ) +
    ggplot2::annotate("text",
                      x = max(X$log2FoldChange) * 0.9, # Adjust x position for label
                      y = -log10(fdr_pvalue_threshold) + 0.1, # Adjust y position for label
                      label = paste0("FDR Adjusted p=", alpha),
                      color = "black",
                      size = 3)
  }

  # Add a label for the raw p-value line
  p <- p + ggplot2::annotate("text",
                            x = max(X$log2FoldChange) * 0.9, # Adjust x position for label
                            y = -log10(alpha) + 0.1, # Adjust y position for label
                            label = paste0("Raw p=", alpha),
                            color = "black",
                            size = 3)

  return(p)
}
