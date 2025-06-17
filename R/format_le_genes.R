#Input: bpn@leading@count_leading_up, bpn@leading@count_leading_down, global_state$data list element
format_le_genes <- function(X, data) {
  X %>%
    rename(Symbol = 1, N = 2) %>%
    arrange(desc(N)) %>%
    inner_join(global_state$map, by = "Symbol") %>%
    inner_join(data, by = "Symbol") %>%
    relocate(Symbol, Name) %>%
    mutate(across(.cols = c(log2FoldChange, pvalue), .fns = ~ round(.x, 3)))
}