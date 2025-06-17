#Input: global state
pca_plot <- function(df) {
  # Extract data components
  counts <- df$counts
  design <- df$design
  
  # Convert counts: set "Symbol" column as row names and transpose the matrix
  pca_matrix <- counts %>%
    tibble::column_to_rownames("Symbol") %>%
    t()
  
  # Perform PCA with centering and scaling
  pca_result <- stats::prcomp(pca_matrix, center = TRUE, scale. = TRUE)
  
  # Generate PCA plot with factoextra; use dplyr to extract the Group variable
  plot <- factoextra::fviz_pca_ind(
    pca_result,
    repel = TRUE,
    habillage = design %>% dplyr::pull(Group),
    addEllipses = TRUE
  ) +
    ggplot2::ggtitle("") +
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
      legend.key.spacing.x = grid::unit(0, "mm"),
      legend.key.spacing.y = grid::unit(0, "mm"),
      legend.box.margin = ggplot2::margin(),
      legend.box = "vertical",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 2, 0, 0), "mm"),
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
  
  label_layer_index <- which(sapply(plot$layers, function(layer)
    inherits(layer$geom, "GeomTextRepel"))) # Find the index of the geom text repel used by factoviz
  
  plot$layers[[label_layer_index]]$aes_params$size <- 9 / ggplot2::.pt # Change the size of the text repel labels
  
  return(plot)
  
}