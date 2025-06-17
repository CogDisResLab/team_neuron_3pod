#Input: table where 1st column is gene symbols to correct
fix_hgnc <- function(X, map) {
  # Validate and correct gene symbols
  corrected <- suppressWarnings(HGNChelper::checkGeneSymbols(dplyr::pull(X, 1), map = map)) %>%
    dplyr::select(Symbol = 1, Suggested.Symbol) %>%
    dplyr::mutate(
      Suggested.Symbol = dplyr::if_else(
        Suggested.Symbol == "" | is.na(Suggested.Symbol),
        Symbol,
        Suggested.Symbol
      ),
      Suggested.Symbol = dplyr::if_else(
        stringr::str_detect(Suggested.Symbol, "///"),
        stringr::word(Suggested.Symbol, 1, sep = "///"),
        Suggested.Symbol
      )
    )
  
  # Join the corrected symbols with the original data and update the gene symbol column
  result <- dplyr::inner_join(corrected, dplyr::rename(X, Symbol = 1), by = "Symbol") %>%
    dplyr::select(-Symbol, Symbol = Suggested.Symbol) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  return(result)
}