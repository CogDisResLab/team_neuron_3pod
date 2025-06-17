#Input: Annotated DrugFindr investigate signature results
format_gene_results <- function(df) {
  df %>%
    dplyr::rename(N = n) %>%
    dplyr::mutate(dplyr::across(.cols = c(log2FoldChange, pvalue), .fns = ~ round(.x, 3)))
}