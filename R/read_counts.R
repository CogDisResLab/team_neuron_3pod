read_counts <- function(file) {
  file %>%
    read_csv() %>%
    rename(Symbol = 1) %>%
    filter(if_all(everything(), ~ . != '')) %>%
    drop_na() %>%
    mutate(Symbol = str_trim(Symbol)) %>%
    mutate(Symbol = if_else(str_detect(Symbol, "///"), str_split_fixed(Symbol, "///", 2)[, 1], Symbol)) %>%
    distinct(Symbol, .keep_all = TRUE)
}