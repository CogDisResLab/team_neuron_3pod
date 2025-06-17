group_count_heatmap <- function(df, useFDR=TRUE, upper_quantile = 0.90, lower_quantile = 0.10, alpha = 0.05) {
  
  # Extract the data element and adjust p-values if requested
  data <- df$data
  
  if (useFDR) {
    data$pvalue <- p.adjust(data$pvalue, method = "fdr")
  }
  
  # Calculate quantiles for log2FoldChange (remove NA if present)
  lower_val <- stats::quantile(data$log2FoldChange, lower_quantile, na.rm = TRUE)
  upper_val <- stats::quantile(data$log2FoldChange, upper_quantile, na.rm = TRUE)
  
  # Filter significant genes based on p-value and fold-change thresholds
  sig_genes <- data %>% 
    dplyr::filter(pvalue <= alpha, 
                  log2FoldChange >= upper_val | log2FoldChange <= lower_val) %>%
    dplyr::pull(Symbol)
  
  if (length(sig_genes) == 0) {
    stop("No significant genes found with the given thresholds when plotting group gene expression heatmap.")
  }
  
  # Get design information for the groups of interest
  groups_of_interest <- unique(c(df$group1, df$group2))
  design_info <- global_state$design %>%
    dplyr::filter(Group %in% groups_of_interest)
  
  # Subset count data for significant genes and samples in the design
  count_matrix <- global_state$counts %>%
    dplyr::select(Symbol, dplyr::all_of(design_info$Sample)) %>%
    dplyr::filter(Symbol %in% sig_genes) %>%
    tibble::column_to_rownames("Symbol") %>%
    as.matrix()
  
  # Scale the matrix by rows
  scaled_mat <- t(scale(t(count_matrix)))
  
  # Determine color scale limits from the scaled matrix
  min_val <- min(scaled_mat, na.rm = TRUE)
  max_val <- max(scaled_mat, na.rm = TRUE)
  
  # Create the heatmap legend with a blue-white-red color ramp
  legend <- ComplexHeatmap::Legend(
    title = NULL,
    col_fun = circlize::colorRamp2(c(min_val, 0, max_val), c("blue", "white", "red")),
    at = c(min_val, 0, max_val),
    labels = c(round(min_val, 2), 0, round(max_val, 2)),
    direction = "vertical",
    labels_gp = grid::gpar(fontsize = 8)
  )
  
  # Create a factor for group annotation and generate a palette with distinct colors
  group_factor <- factor(design_info$Group)
  group_palette <- randomcoloR::distinctColorPalette(nlevels(group_factor))
  names(group_palette) <- levels(group_factor)
  
  # Set up heatmap annotation for group information
  annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_factor,
    show_annotation_name = FALSE,
    col = list(Group = group_palette),
    annotation_legend_param = list(
      Group = list(
        color_bar = "discrete",
        title = NULL,
        direction = "horizontal",
        title_position = "topcenter",
        title_gp = grid::gpar(fontsize = 8),
        labels_gp = grid::gpar(fontsize = 8),
        nrow = 1
      )
    )
  )
  
  # Create the heatmap without row names or dendrogram, rotating column names for clarity
  heatmap_plot <- ComplexHeatmap::Heatmap(
    scaled_mat,
    top_annotation = annotation,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = FALSE
  )
  
  # Draw the heatmap with specified padding and legend settings
  ComplexHeatmap::draw(
    heatmap_plot,
    padding = grid::unit(c(2, 15, 2, 2), "mm"),
    heatmap_legend_list = legend,
    heatmap_legend_side = "right",
    annotation_legend_side = "top"
  )
  
}