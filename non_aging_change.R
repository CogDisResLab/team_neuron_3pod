# Idenfity the genes that have a non-aging change

library(tidyverse)

deg_files <- list.files("extdata", "clean") |>
    set_names(~ str_remove(.x, ".csv"))


filter_pngf <- function(file_a, file_b) {
    data_a <- read_csv(file.path("extdata", file_a)) |>
        rename(LFC_pngf = logFC, PVal_pngf = P.Value)

    data_b <- read_csv(file.path("extdata", file_b)) |>
        rename(LFC_wt = logFC, PVal_wt = P.Value)

    data_combined <- left_join(data_a, data_b, by = "gene_name")
}


pairs <- tribble(
    ~file_a, ~file_b,
    deg_files[1L], deg_files[4L],
    deg_files[2L], deg_files[5L],
    deg_files[3L], deg_files[6L]
)

comparisons <- pmap(pairs, filter_pngf) |>
    map2(
        c("1m_12m", "1m_3m", "3m_12m"),
        ~ write_csv(.x, file.path("results", str_glue("comparison_{.y}.csv")))
    )
