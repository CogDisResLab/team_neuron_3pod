# Idenfity the genes that have a non-aging change

library(tidyverse)

deg_files <- list.files("extdata", "clean") |>
  set_names(~ str_remove(.x, ".csv"))


filter_pngf <- function(file_a, file_b) {
  data_a <- read_csv(file.path("extdata", file_a)) |>
    rename(LFC_pngf = logFC, PVal_pngf = P.Value) |>
    filter(PVal_pngf < 0.05)

  data_b <- read_csv(file.path("extdata", file_b)) |>
    rename(LFC_wt = logFC, PVal_wt = P.Value) |>
    filter(PVal_wt < 0.05)

  data_combined <- data_a |>
    left_join(data_b, by = "gene_name")
}


pairs <- tribble(
  ~file_a, ~file_b,
  deg_files[1], deg_files[4],
  deg_files[2], deg_files[5],
  deg_files[3], deg_files[6]
)

comparisons <- pmap(pairs, filter_pngf) |>
  map2(
    c("1m_12m", "1m_3m", "3m_12m"),
    ~ write_csv(.x, file.path("results", str_glue("comparison_{.y}.csv")))
  )
