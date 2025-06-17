#' Read and Clean Differential Expression Gene Results
#'
#' This function reads a CSV file generated from differential expression testing analyses and cleans the data for downstream use.
#' The input CSV is expected to have at least three columns:
#' \describe{
#'   \item{Column 1}{Standardized gene symbols (e.g., HGNC, RGD, MGI).}
#'   \item{Column 2}{Log2 fold change values (numeric, with both positive and negative values).}
#'   \item{Column 3}{Raw p-values (numeric values between 0 and 1).}
#' }
#'
#' The function performs several sanity checks including:
#' \itemize{
#'   \item Verifying file existence and proper column count.
#'   \item Ensuring that the gene symbols are character strings, log2 fold changes are numeric with a mix of positive and negative values,
#'         and p-values are numeric within the [0, 1] interval.
#'   \item Checking for a high duplicate rate in the gene symbols.
#'   \item Cleaning the data by removing rows with blanks or NAs, squishing extra whitespace, and, if multiple symbols are present,
#'         retaining only the first.
#'   \item Comparing row counts before and after cleaning.
#'   \item Additional quality checks on the cleaned data including:
#'         \item Filtering for significant features (p-value <= 0.05) and issuing a warning if fewer than 100 rows remain.
#'         \item Warning if a large majority (over 75%) of log2 fold change values lie between -0.25 and 0.25.
#'         \item Detecting infinite or extreme log2 fold change values.
#'         \item Evaluating the p-value distribution and overall variance in log2 fold changes.
#'         \item Warning if the final cleaned dataset has too few rows.
#' }
#'
#' @param file A character string specifying the path to the CSV file containing the differential expression results.
#'
#' @return A tibble (data frame) with three columns: \code{Symbol}, \code{log2FoldChange}, and \code{pvalue}. 
#'         Warnings are issued if potential formatting or quality issues are detected.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'deg_results.csv' is a CSV file from a differential expression analysis:
#'   result_df <- read_deg("deg_results.csv")
#' }
#'
#' @importFrom readr read_csv
#' @importFrom dplyr select rename filter mutate distinct n_distinct if_all
#' @importFrom tidyr drop_na
#' @importFrom stringr str_squish str_detect word
#'
#' @export
read_deg <- function(file) {
  #### 1. INITIAL FILE & COLUMN VALIDATION ####
  
  # Check if the file exists
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  # Attempt to read the file, capturing any errors
  df <- tryCatch(
    {
      readr::read_csv(file, 
                      show_col_types = FALSE,
                      name_repair = "unique_quiet")
    },
    error = function(e) {
      stop("Error reading file: ", e$message)
    }
  )
  
  # Verify the data has at least three columns
  if (ncol(df) < 3) {
    stop("Input file must contain at least 3 columns.")
  }
  
  # Select and rename the first three columns
  df <- df %>%
    dplyr::select(1:3) %>%
    dplyr::rename(Symbol = 1, log2FoldChange = 2, pvalue = 3)
  
  # Basic sanity checks on data types and expected ranges:
  
  # Check that 'Symbol' is a character vector.
  if (!is.character(df$Symbol)) {
    stop("Column 'Symbol' must be of type character.")
  }
  
  # Check that 'log2FoldChange' is numeric and has both positive and negative values.
  if (!is.numeric(df$log2FoldChange)) {
    stop("Column 'log2FoldChange' must be numeric.")
  }
  if (min(df$log2FoldChange, na.rm = TRUE) >= 0 || max(df$log2FoldChange, na.rm = TRUE) <= 0) {
    stop("Column 'log2FoldChange' must contain both positive and negative values.")
  }
  
  # Check that 'pvalue' is numeric and all values lie between 0 and 1.
  if (!is.numeric(df$pvalue)) {
    stop("Column 'pvalue' must be numeric.")
  }
  if (min(df$pvalue, na.rm = TRUE) < 0 || max(df$pvalue, na.rm = TRUE) > 1) {
    stop("Column 'pvalue' values must range between 0 and 1.")
  }
  
  #### 2. PRE-CLEANING CHECKS ####
  
  # Check for duplicate gene symbols in raw data (ignoring blanks and NAs)
  raw_df <- df %>% dplyr::filter(Symbol != "", !is.na(Symbol))
  dup_fraction <- 1 - (dplyr::n_distinct(raw_df$Symbol) / nrow(raw_df))
  if (dup_fraction > 0.15) {
    warning("High duplicate rate: ", round(dup_fraction * 100, 2), 
            "% of gene symbols are duplicates in the raw data.")
  }
  
  # Store the original row count before cleaning
  original_rows <- nrow(df)
  
  #### 3. DATA CLEANING & POST-CLEANING CHECKS ####
  
  # Clean and transform data:
  df <- df %>%
    # Remove rows where any column is blank and drop NA values
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != '')) %>%
    tidyr::drop_na() %>%
    # Clean up 'Symbol': squish whitespace and, if multiple symbols exist (separated by "///"),
    # retain only the first one.
    dplyr::mutate(
      Symbol = stringr::str_squish(Symbol),
      Symbol = dplyr::if_else(
        stringr::str_detect(Symbol, "///"),
        stringr::word(Symbol, 1, sep = "///"),
        Symbol
      ),
      Symbol = stringr::str_squish(Symbol)
    ) %>%
    # Remove duplicate symbols
    dplyr::distinct(Symbol, .keep_all = TRUE)
  
  # Compare row count before and after cleaning.
  final_rows <- nrow(df)
  if (final_rows != original_rows) {
    warning("The number of rows has changed from ", original_rows, " to ", final_rows, 
            ". This may indicate duplicates, NAs, or blank values were present in the data.")
  }
  
  # Additional sanity checks on the cleaned data:
  
  ## (a) Significant Features Check: Warn if few rows remain after filtering for pvalue <= 0.05.
  significant_df <- df %>% dplyr::filter(pvalue <= 0.05)
  if (nrow(significant_df) < 100) {
    warning("After filtering for pvalue <= 0.05, the number of rows is ", nrow(significant_df), 
            ", which is less than 100. This may indicate that few features are statistically significant.")
  }
  
  ## (b) Minimal Change Check: Warn if a large majority of log2FoldChange values are near zero.
  num_in_range <- sum(df$log2FoldChange >= -0.25 & df$log2FoldChange <= 0.25, na.rm = TRUE)
  proportion_in_range <- num_in_range / nrow(df)
  if (proportion_in_range > 0.75) {
    warning("More than 75% of log2FoldChange values are between -0.25 and 0.25, implying few differences.")
  }
  
  ## (c) Extreme Values Check for log2FoldChange:
  if (any(is.infinite(df$log2FoldChange))) {
    warning("Infinite values found in the log2FoldChange column.")
  }
  extreme_threshold <- 10
  prop_extreme <- sum(abs(df$log2FoldChange) > extreme_threshold, na.rm = TRUE) / nrow(df)
  if (prop_extreme > 0.10) {
    warning("More than 10% of log2FoldChange values exceed an absolute value of ", extreme_threshold, 
            ". Please check for possible calculation errors or formatting issues.")
  }
  
  ## (d) P-value Distribution Check: Warn if an unusually high proportion of pvalues are exactly 0 or 1.
  prop_p0 <- sum(df$pvalue == 0, na.rm = TRUE) / nrow(df)
  prop_p1 <- sum(df$pvalue == 1, na.rm = TRUE) / nrow(df)
  if (prop_p0 > 0.5) {
    warning("More than 50% of p-values are 0. This may indicate an issue with the statistical tests or data formatting.")
  }
  if (prop_p1 > 0.5) {
    warning("More than 50% of p-values are 1. This may indicate an issue with the statistical tests or data formatting.")
  }
  
  ## (e) Low Variance Check for log2FoldChange: Warn if the variance is very low.
  if (var(df$log2FoldChange, na.rm = TRUE) < 0.01) {
    warning("Variance of log2FoldChange is very low, suggesting there may be minimal differences between conditions.")
  }
  
  ## (f) Near-Empty Data Check: Warn if the cleaned data contains too few rows.
  if (final_rows < 10) {
    warning("The cleaned data contains fewer than 10 rows. Please check the input file for formatting issues.")
  }
  
  return(df)
}