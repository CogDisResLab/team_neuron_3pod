#Input: Vector of Gene symbols
#Output: Tibble of enrichR results
do_enrichr <- function(symbols, alpha = 0.05) {
  
  Sys.sleep(5) # Enforce a 5 second pause to avoid known issues with enrichR api backend returning duplicates if you request too fast
  
  # TODO: re-implement enrichR locally using the GMT file? 
  
  dbs = c("GO_Biological_Process_2025", 
          "GO_Molecular_Function_2025", 
          "GO_Cellular_Component_2025")
  
  columns = c("Biological_Process",
              "Molecular_Function",
              "Cellular_Component")
  
  res <- quiet(symbols %>%
          enrichR::enrichr(databases = dbs) %>% 
          purrr::map2_dfr(columns, ~ mutate(.x, namespace = .y)) %>%
          dplyr::filter(Adjusted.P.value <= alpha) %>%
          tidyr::extract(Term, "GOID", "(GO:\\d+)", remove = FALSE) %>%
          tidyr::extract(Term, "Term", "(.*?)\\(GO:\\d+\\)"))
  res

}