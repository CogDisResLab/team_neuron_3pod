# Extract and dump lincs perturbagen data

library(tidyverse)
library(readxl)

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

    combined_drugs <- bind_rows(
        data_to_process[[sheets[[1L]]]],
        data_to_process[[sheets[[2L]]]]
    )
    write_csv(
        combined_drugs,
        file.path(result_path, str_glue("lincs_{comp_name}_drugs.csv"))
    )
}

process_names <- function(name) {
    genotype <- str_detect(name, fixed("ProNGF"))
    time_point <- str_detect(name, fixed("X1"))

    case_when(
        genotype & time_point ~ "PNGF_1M",
        genotype & !time_point ~ "PNGF_3M",
        !genotype & time_point ~ "WT_1M",
        !genotype & !time_point ~ "WT_3M",
    )
}

process_dataset <- function(filename, comparison) {
    filepath <- file.path("results", "lincs", filename)

    read_csv(filepath) |>
        select(starts_with("Target"), FDA_Phase = max_fda_phase, Similarity) |>
        mutate({{ comparison }} := Similarity) |>
        select(-Similarity)
}

load(file.path("results", "report_environment.RData"))

lincs_results <- global_state$data |>
    imap(~ extract_lincs_data(.y, .x))


combined_results <- list.files(file.path("results", "lincs"), "csv", recursive = TRUE) |>
    keep(~ str_detect(.x, fixed("drugs"))) |>
    set_names(~ process_names(.x)) |>
    imap(process_dataset) |>
    reduce(full_join) |>
    write_csv(file.path("results", "combined_lincs_similarity.csv"))
