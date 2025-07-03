# Load libraries
library(tidyverse)

# Define genes of interest
genes_of_interest <- c("Pion", "Kcnab2", "Neurod1", "Neurod2", "Gnai")

# Discover files and name them by comparison
deg_files <- list.files("extdata", pattern = "clean.*\\.csv$") |>
  discard(~ str_detect(.x, "1m-12m")) |>
  set_names(~ str_extract(.x, "^[^\\-]+"))

# Function to process individual DEG files
process_data <- function(filename) {
  filepath <- file.path("extdata", filename)
  read_csv(filepath, show_col_types = FALSE) |>
    mutate(
      Adj.P.Value = p.adjust(P.Value, method = "fdr"),
      logFC = round(logFC, 6)
    ) |>
    select(gene_name, logFC, Adj.P.Value) |>
    filter(gene_name %in% genes_of_interest)
}

# Combine and annotate data
combined_data <- map(deg_files, process_data) |>
  bind_rows(.id = "Comparison") |>
  mutate(
    Group = case_when(
      str_detect(Comparison, "1m") ~ "1M3M",
      str_detect(Comparison, "3m") ~ "3M12M"
    ),
    Genotype = case_when(
      str_starts(Comparison, "wt") ~ "Wild Type",
      str_starts(Comparison, "pngf") ~ "proNGF"
    ),
    Comparison = factor(Comparison, levels = c("wt_1m", "wt_3m", "pngf_1m", "pngf_3m")),
    gene_name = factor(gene_name, levels = rev(genes_of_interest)),
    sig_label = if_else(Adj.P.Value < 0.05, "*", "")
  )

# Two-color palette
genotype_colors <- c(
  "Wild Type" = "#1b9e77",
  "proNGF"    = "#d95f02"
)

# Final plot
ggplot(combined_data, aes(x = logFC, y = gene_name, fill = Genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.75, aes(group = Comparison)) +
  geom_text(
    aes(label = sig_label),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 5,
    color = "black"
  ) +
  facet_wrap(~Group, nrow = 1, labeller = as_labeller(c("1M3M" = "1 to 3 Months", "3M12M" = "3 to 12 Months"))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
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
  theme_bw(base_size = 16) + # sets base size for all text elements
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggsave("lfc_pyramid_plot.png", width = 11, height = 8.5, dpi = 300, path = "figures")
