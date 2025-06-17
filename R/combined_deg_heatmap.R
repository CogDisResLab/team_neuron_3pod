#Input: a list of DEG tables
combined_deg_heatmap <- function(df,
                                 num_genes = 50,
                                 useFDR = TRUE,
                                 alpha = 0.05) {
  if (useFDR) {
    df <- df %>%
      purrr::map(~ dplyr::mutate(., pvalue = p.adjust(pvalue, method = "fdr")))
  }
  
  # Identify top DEGs in each dataset
  top_degs <- df %>%
    purrr::map(
      ~ dplyr::filter(., pvalue <= alpha) %>%
        dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
        dplyr::slice_head(n = num_genes) %>%
        dplyr::pull(Symbol)
    ) %>%
    purrr::flatten_chr() %>%
    unique()
  
  mat <- df %>%
    purrr::map( ~ dplyr::filter(., Symbol %in% top_degs) %>%
                  dplyr::select(Symbol, log2FoldChange)) %>%
    dplyr::bind_rows(.id = "Group") %>%
    tidyr::pivot_wider(names_from = Group, values_from = log2FoldChange) %>%
    tidyr::drop_na() %>%
    tibble::column_to_rownames(var = "Symbol") %>%
    as.matrix()
  
  min = min(mat, na.rm = T)
  max = max(mat, na.rm = T)
  
  lgd = ComplexHeatmap::Legend(
    title = "log2FC",
    col_fun = circlize::colorRamp2(c(min, 0, max), c("blue", "white", "red")),
    at = c(min, 0, max),
    labels = c(round(min, 2), 0, round(max, 2)),
    direction = "vertical",
    labels_gp = grid::gpar(fontsize = 7),
    title_gp = grid::gpar(fontsize = 7)
  )
  
  ht = ComplexHeatmap::Heatmap(
    mat,
    column_names_rot = 0,
    column_names_gp = grid::gpar(fontsize = 9),
    column_names_centered = TRUE,
    row_names_gp = grid::gpar(fontsize = 5),
    show_row_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE
  )
  
  ComplexHeatmap::draw(
    ht,
    padding = grid::unit(c(0, 0, 0, 0), "mm"),
    heatmap_legend_list = lgd,
    heatmap_legend_side = "right"
  )
  
}