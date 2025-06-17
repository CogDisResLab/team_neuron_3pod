#Input: Annotated DrugFindr investigate signature results
generate_moa_report <- function(annotated_signature) {
  moa_report <- annotated_signature %>%
    dplyr::select(integratedMoas, Target, GeneTargets) %>%
    dplyr::filter(integratedMoas != "" & !is.na(integratedMoas)) %>%
    tidyr::separate_rows(integratedMoas, sep = "\\|") %>%
    tidyr::separate_rows(GeneTargets, sep = "\\|") %>%
    dplyr::mutate(dplyr::across(c(integratedMoas, GeneTargets), stringr::str_trim)) %>%
    dplyr::group_by(integratedMoas) %>%
    dplyr::summarise(
      Target = paste(unique(Target), collapse = "|"),
      GeneTargets = paste(unique(GeneTargets), collapse = "|"),
      N = dplyr::n()
    ) %>%
    dplyr::arrange(dplyr::desc(N))
    
    moa_report
}
