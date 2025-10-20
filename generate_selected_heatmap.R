# Generate a specific heatmap with unnecessary pathways removed

suppressPackageStartupMessages({
    library(PAVER)
    library(tidyverse)
})

load(file.path("results", "report_environment.RData"))

paver_results <- global_state$results

enrichr_paver <- paver_results$ENRICHR_PAVER_result

selected_themes <- c(
    "regulation of cytoplasmic translation",
    "calcium ion transmembrane transporter activity",
    "regulation of exocytosis", "nervous system development",
    "mitochondrial ATP synthesis coupled electron transport",
    "asymmetric synapse"
)

modified_clusters <- filter(
    enrichr_paver[["clustering"]],
    Cluster %in% selected_themes
)

enrichr_paver$clustering <- modified_clusters

paver_plot <- PAVER_hunter_plot(enrichr_paver)


ggsave("selected_pathway_heatmap.png",
    bg = "white",
    width = 7.5 * 3L,
    height = 10L * 3L,
    path = "figures"
)

ggsave("selected_pathway_heatmap.svg",
    width = 7.5 * 3L,
    height = 10L * 3L,
    path = "figures"
)
