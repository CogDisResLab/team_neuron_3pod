# Check the correlation between the DEGs

library(tidyverse)

deg_files <- list.files("extdata", "clean") |>
  set_names(~ str_remove(.x, ".csv"))

perform_correlation <- function(file_a, file_b) {
  data_a <- read_csv(file.path("extdata", file_a))
  data_b <- read_csv(file.path("extdata", file_b))

  correlation <- cor(data_a$logFC, data_b$logFC)
}

pairs <- tribble(
  ~file_a, ~file_b,
  deg_files[1], deg_files[4],
  deg_files[2], deg_files[5],
  deg_files[3], deg_files[6]
)

cors <- pmap(pairs, perform_correlation)
