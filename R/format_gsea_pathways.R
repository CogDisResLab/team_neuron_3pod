#Input: bpn@gsea_sig_pos_enriched, bpn@gsea_sig_neg_enriched
format_gsea_pathways <- function(df) {
  df %>%
    dplyr::mutate(
      Term = stringr::str_remove(pathway, "%.*") %>% stringr::str_to_sentence(),
      padj = round(padj, 3),
      NES = round(NES, 3), 
      leadingEdge, 
      .keep = "none") %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    dplyr::relocate(Term)
}