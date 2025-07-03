# Load libraries
library(tidyverse)

# Define genes of interest
genes_of_interest <- c("Pion", "Kcnab2", "Neurod1", "Neurod2", "Gnai")

# Discover files and name them by comparison
deg_files <- list.files("extdata", pattern = "clean.*\\.csv$") |>
    discard(~ str_detect(.x, fixed("1m-12m"))) |>
    set_names(~ str_extract(.x, "^[^\\-]+")) # nolint: nonportable_path_linter.

# Function to process individual DEG files
process_data <- function(filename) {
    filepath <- file.path("extdata", filename)
    read_csv(filepath, show_col_types = FALSE) |>
        mutate(
            Adj.P.Value = p.adjust(P.Value, method = "fdr"),
            logFC = round(logFC, 6L)
        ) |>
        select(gene_name, logFC, Adj.P.Value) |>
        filter(gene_name %in% genes_of_interest)
}

# Combine and annotate data
combined_data <- map(deg_files, process_data) |>
    bind_rows(.id = "Comparison") |>
    mutate(
        Group = case_when(
            str_detect(Comparison, fixed("1m")) ~ "1M3M",
            str_detect(Comparison, fixed("3m")) ~ "3M12M"
        ),
        Genotype = case_when(
            str_starts(Comparison, fixed("wt")) ~ "Wild Type",
            str_starts(Comparison, fixed("pngf")) ~ "proNGF"
        ),
        Comparison = factor(Comparison, levels = c("wt_1m", "wt_3m", "pngf_1m", "pngf_3m")),
        gene_name = factor(gene_name, levels = rev(genes_of_interest)),
        sig_label = if_else(Adj.P.Value < 0.05, "*", "")
    )

# Two-color palette
genotype_colors <- c(
    "Wild Type" = "#1b9e77",
    proNGF = "#d95f02"
)

# Final plot
ggplot(combined_data, aes(x = logFC, y = gene_name, fill = Genotype)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.75, aes(group = Comparison)) +
    geom_text(
        aes(label = sig_label),
        position = position_dodge(width = 0.8),
        hjust = -0.2,
        size = 5L,
        color = "black"
    ) +
    facet_wrap(~Group, nrow = 1L, labeller = as_labeller(c("1M3M" = "1 to 3 Months", "3M12M" = "3 to 12 Months"))) +
    geom_vline(xintercept = 0L, color = "black", linetype = "dashed") +
    scale_fill_manual(values = genotype_colors) +
    scale_x_continuous(
        name = "log₂ Fold Change",
        breaks = seq(-1.5, 1.5, 0.5),
        limits = c(-1.5, 1.5)
    ) +
    ylab("") +
    labs(
        fill = "Genotype",
        title = "Differential Gene Expression Across Age and Genotype"
    ) +
    theme_bw(base_size = 16L) + # sets base size for all text elements
    theme(
        legend.position = "bottom",
        legend.title = element_text(size = 16L, face = "bold"),
        legend.text = element_text(size = 14L),
        axis.title.x = element_text(size = 18L, face = "bold"),
        axis.text.x = element_text(size = 14L),
        axis.text.y = element_text(size = 14L),
        strip.text = element_text(size = 16L, face = "bold"),
        plot.title = element_text(size = 20L, face = "bold", hjust = 0.5)
    )

ggsave("lfc_pyramid_plot.png", width = 11L, height = 8.5, dpi = 300L, path = "figures")
