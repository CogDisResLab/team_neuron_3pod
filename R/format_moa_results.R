#Input: Annotated DrugFindr investigate signature results
format_moa_results <- function(df) {
  df %>%
    dplyr::select(MOA = integratedMoas, N)
  
}