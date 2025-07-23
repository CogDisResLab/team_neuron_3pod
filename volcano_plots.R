# Generate volcano plots that are more figure-ready

library(tidyverse)

deg_files <- list.files("extdata", "csv") |>
    keep(~ str_detect(.x, fixed("clean"))) |>
    discard(~ str_detect(.x, fixed("1m-12m"))) |>
    set_names(~ str_extract(.x, "(^\\w+)-", 1L)) # nolint: nonportable_path_linter.

process_data <- function(filename) {
    filepath <- file.path("extdata", filename)
    deg <- read_csv(filepath) |>
        mutate(
            Adj.P.Value = stats::p.adjust(P.Value, method = "fdr"),
            Significant = Adj.P.Value <= 0.05,
            Regulation = sign(logFC),
            Color = case_when(
                Significant & Regulation == 1L ~ "darkred",
                Significant & Regulation == -1L ~ "darkblue",
                .default = "black"
            ),
            X = logFC,
            Y = -log10(Adj.P.Value),
            Score = abs(X) * Y,
        ) |>
        arrange(desc(Score)) |>
        select(gene_name, X, Y, Color, Significant)
}

genes_of_interest <- c()

combined_data <- map(deg_files, process_data) |>
    bind_rows(.id = "Comparison") |>
    mutate(
        Facet = case_when(
            Comparison == "pngf_1m" ~ "proNGF 1m vs 3m",
            Comparison == "pngf_3m" ~ "proNGF 3m vs 12m",
            Comparison == "wt_1m" ~ "Wild Type 1m vs 3m",
            Comparison == "wt_3m" ~ "Wild Type 3m vs 12m",
        ),
        Label = if_else(gene_name %in% genes_of_interest, gene_name, NA_character_)
    ) |>
    select(-Comparison, -gene_name)

g <- ggplot(combined_data, aes(
    x = X, y = Y,
    color = Color, label = Label
))

threshold <- -log10(0.05)

p <- g +
    geom_hline(yintercept = threshold, lty = "dashed", lwd = 0.25) +
    geom_vline(xintercept = c(-1L, 1L), lwd = 0.5, color = "grey") +
    geom_point() +
    scale_color_manual(
        limits = c("darkred", "black", "darkblue"),
        labels = c("Upregulated", "NS", "Downregulated"),
        values = c("darkred", "black", "darkblue"),
        name = ""
    ) +
    ggrepel::geom_text_repel(
        max.overlaps = Inf,
        segment.color = "grey30",
        box.padding = 1L,
        segment.size = 2L,
        min.segment.length = 1L,
        show.legend = FALSE
    ) +
    scale_x_continuous(name = expression(log["2"]("FC")), breaks = seq(-10L, 10L, 1L)) +
    scale_y_continuous(name = expression(-log["10"](italic(p) * "-value")), breaks = seq(0L, 10L, 1L)) +
    facet_wrap(facets = "Facet", nrow = 2L, ncol = 2L, axes = "all", dir = "v") +
    theme_minimal(base_size = 36L) +
    guides(colour = guide_legend(override.aes = list(size = 12L))) +
    theme(legend.position = "bottom", legend.box.margin = margin(t = 10L))

p

ggsave("combined_volcano_plots.png",
    bg = "white",
    height = 7.5 * 3L,
    width = 10L * 3L,
    path = "figures"
)

ggsave("combined_volcano_plots.svg",
    bg = "white",
    height = 7.5 * 3L,
    width = 10L * 3L,
    path = "figures"
)
