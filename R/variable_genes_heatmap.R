#Input: global_state
variable_genes_heatmap <- function(df, num_genes = 500) {
  # Convert counts data frame to a matrix with gene symbols as rownames
  mat <- df$counts %>%
    tibble::column_to_rownames(var = "Symbol") %>%
    as.matrix()
  
  # Calculate the variance of each gene (row) and order genes by decreasing variance
  indices <- apply(mat, 1, var) %>% order(decreasing = TRUE)
  
  # Ensure we do not request more genes than available
  num_genes <- min(num_genes, nrow(mat))
  
  # Subset to the top 'num_genes' and scale each row
  mat <- mat[indices[1:num_genes], , drop = FALSE]
  mat <- t(scale(t(mat)))
  
  # Determine the minimum and maximum for color scaling
  min_val <- min(mat, na.rm = TRUE)
  max_val <- max(mat, na.rm = TRUE)
  
  # Create a legend for the heatmap
  lgd <- ComplexHeatmap::Legend(
    title = NULL,
    col_fun = circlize::colorRamp2(c(min_val, 0, max_val), c("blue", "white", "red")),
    at = c(min_val, 0, max_val),
    labels = c(round(min_val, 2), 0, round(max_val, 2)),
    direction = "vertical",
    labels_gp = grid::gpar(fontsize = 8)
  )
  
  # Process group annotation data
  group_data <- as.factor(df$design$Group)
  group_cols <- randomcoloR::distinctColorPalette(k = nlevels(group_data)) %>%
    purrr::set_names(levels(group_data))
  
  annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_data,
    show_annotation_name = FALSE,
    col = list(Group = group_cols),
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
  
  # Create and draw the heatmap
  ht <- ComplexHeatmap::Heatmap(
    mat,
    top_annotation = annotation,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = FALSE
  )
  
  ComplexHeatmap::draw(
    ht,
    padding = grid::unit(c(2, 15, 2, 2), "mm"),
    heatmap_legend_list = lgd,
    heatmap_legend_side = "right",
    annotation_legend_side = "top"
  )
}