#Input: bpn@enrichr@sig_up_enrichr, bpn@enrichr@sig_down_enrichr
format_enrichr_pathways <- function(df) {
  df %>%
    mutate(
      Term = stringr::str_remove(Term, "\\(([^)]+)\\)") %>% stringr::str_to_sentence(),
      padj = round(Adjusted.P.value, 3),
      CS = round(Combined.Score, 3),
      Genes,
      .keep = "none"
    ) %>%
    dplyr::arrange(dplyr::desc(CS), padj) %>%
    dplyr::relocate(padj, CS, .after = Term)
}