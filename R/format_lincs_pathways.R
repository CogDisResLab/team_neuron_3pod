#Input: .$results$lincs$concordant_pathways, .$results$lincs$discordant_pathways
format_lincs_pathways <- function(X) {
  X %>%
    arrange(desc(Combined.Score), Adjusted.P.value) %>%
    select(Term, padj = Adjusted.P.value, CS = Combined.Score, Genes) %>%
    mutate(
      Term = gsub("\\(([^)]+)\\)", "", Term) %>% str_to_sentence,
      padj = round(padj, 3),
      CS = round(CS, 3)
    )
}