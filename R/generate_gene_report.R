#Input: Annotated DrugFindr investigate signature results | signature data
generate_gene_report <- function(annotated_signature, signature) {
  gene_report <- annotated_signature %>%
    dplyr::select(Symbol = GeneTargets) %>%
    dplyr::filter(Symbol != "" & !is.na(Symbol)) %>%
    tidyr::separate_rows(Symbol, sep = "\\|") %>%
    dplyr::mutate(Symbol = stringr::str_trim(Symbol)) %>%
    dplyr::count(Symbol) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::inner_join(signature, by = "Symbol") %>% 
    dplyr::inner_join(global_state$hgnc, by = "Symbol")
    
  gene_report
}