# Extract and dump lincs perturbagen data

library(tidyverse)

extract_lincs_data <- function(comparison, lincs_data) {
  comp_name <- make.names(comparison)
  data_to_process <- pluck(
    lincs_data,
    "results", "lincs"
  )

  result_path <- file.path("results", "lincs", comp_name)

  dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

  sheets <- c(
    "concordant", "discordant",
    "concordant_moa_report", "discordant_moa_report"
  )

  for (sheet in sheets) {
    filename <- str_glue("lincs_{comp_name}_{sheet}.csv")
    filepath <- file.path(result_path, filename)
    write_csv(data_to_process[[sheet]], filepath)
  }
}

global_state <- read_rds("global_state.RDS")

lincs_results <- global_state$data |>
  imap(~ extract_lincs_data(.y, .x))
