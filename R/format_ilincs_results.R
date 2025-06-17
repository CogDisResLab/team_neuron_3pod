#Input: Annotated DrugFindr investigate signature results
format_ilincs_results <- function(df) {
  df %>%
    dplyr::mutate(Similarity = round(Similarity, 3)) %>%
    dplyr::select(
      Similarity,
      Perturbagen = Target,
      MOA = integratedMoas,
      GeneTargets,
      Cell = TargetCellLine,
      Tissue = tissue,
      Concentration = TargetConcentration,
      Time = TargetTime,
      `FDA Phase` = max_fda_phase)
}  