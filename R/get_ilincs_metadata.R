#Input: Vector of LINCS signature ids from drugFindr, e.g. LINCSCP_136625
#Output: Tibble of iLINCS signature metadata
get_ilincs_metadata <- function(SignatureID) {
  
  url <- "https://www.ilincs.org/api/SignatureMeta/findMany"
  body <- list(signatures = jsonlite::toJSON(SignatureID))
  
  metadata <- httr::POST(url, body = body, encode = "json") %>%
    httr::content(as = "text") %>%
    jsonlite::fromJSON() %>%
    purrr::pluck("data") %>%
    select(TargetSignature = signatureid, 
           tissue, 
           integratedMoas, 
           GeneTargets) #where(~ any(!is.na(.x)))
  
  metadata
}