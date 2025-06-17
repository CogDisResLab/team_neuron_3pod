# Load libraries and set options
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(magrittr)
  
  library(jsonlite)
  library(httr)
  library(knitr)
 
  library(ggupset)
  library(ggrepel)
  library(ggVennDiagram)
  
  library(DT)
  library(circlize)
  library(randomcoloR)
  library(ComplexHeatmap)
  library(factoextra)

  library(HGNChelper)
  library(BioPathNet)
  library(drugfindR)
  library(PAVER)
  library(enrichR)
  library(babelgene)
})

suppressMessages(
  options(
    readr.show_col_types = FALSE,
    timeout = 999,
    rlib_name_repair_verbosity = "quiet"
  )
)

set.seed(123)

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  echo = FALSE,
  out.width = "100%",
  out.height = "100%",
  fig.align = "center",
  dpi = 300 # TODO: 600, or use vector plots?
)

# Load functions used throughout the report
list.files(here::here("R"), full.names = TRUE) %>%
  purrr::walk(source)

# Setup report state
configuration_yml <- here::here("configuration.yml")

if(!file.exists(configuration_yml)) {
  stop("Configuration file not found")
}

global_state <- yaml::read_yaml(file = configuration_yml, readLines.warn = FALSE)

# This is the HGNCHelper precomputed mapping file
global_state$humanmap <- here::here("extdata", "assets", global_state$humanmap) %>%
  readr::read_csv()

#These are gene annotations for each species
global_state$hgnc <- here::here("extdata", "assets", global_state$hgnc) %>%
  readr::read_csv() %>%
  dplyr::select(Symbol = symbol, Name = name)

global_state$mgi <- here::here("extdata", "assets", global_state$mgi) %>%
  readr::read_csv() %>%
  dplyr::select(Symbol = `Marker Symbol`, Name = `Marker Name`)

global_state$rgd <- here::here("extdata", "assets", global_state$rgd) %>%
  readr::read_csv() %>%
  dplyr::select(Symbol = SYMBOL, Name = NAME)

global_state$map <- switch(
  global_state$species,
  "human" = global_state$hgnc,
  "mouse" = global_state$mgi,
  "rat"   = global_state$rgd,
  stop("Unsupported species. Select one of 'human', 'mouse', or 'rat' in the configuration file.")
)

#This is the scraped LINCS FDA data
global_state$lincs_fda <- here::here("extdata", "assets", global_state$lincs_fda) %>%
  read_csv()

#This is the GMT file
global_state$gmt <- here::here("extdata", "assets", global_state$gmt) 

#Theses are files for PAVER
global_state$embeddings <- here::here("extdata", "assets", global_state$embeddings) %>%
  readr::read_rds()

global_state$term2name <- here::here("extdata", "assets", global_state$term2name) %>%
  read_rds()

#Read counts and design if specified
design_file <- here::here("extdata", global_state$design)
counts_file <- here::here("extdata", global_state$counts)

if(!rlang::is_empty(global_state$design) && !rlang::is_empty(global_state$counts) && file.exists(design_file) && file.exists(counts_file)) {
  global_state$design <- readr::read_csv(design_file)
  global_state$counts <- read_counts(counts_file)
  
  if (global_state$species == "human") {
    #Correct HGNC symbols of the counts if human data
    global_state$counts <- fix_hgnc(global_state$counts, global_state$humanmap)
  }
  
  global_state$using_counts <- TRUE
  
} else {
  global_state$using_counts <- FALSE
}

# Load input DEG data
global_state$data <- global_state$data %>%
  purrr::map(function(x) {
   
    x$name  = stringr::str_c(x$group1, " vs ", x$group2)
    x$data  = read_deg(here::here("extdata", x$file))
    x$results = list()
    
    if (global_state$species == "human") {
      x$data <- fix_hgnc(x$data, global_state$humanmap)
    }
    
    x$bpn <- BioPathNet::prepare_data(
      x$data$Symbol, x$data$log2FoldChange, x$data$pvalue
    )
    
    x
  })

global_state$data <- purrr::set_names(global_state$data, purrr::map_chr(global_state$data, "name"))

#Create global results store
global_state$results <- list()