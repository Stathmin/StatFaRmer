# Processed Data Generation
# This generates the processed data files with logit transforms and outlier handling

# cleaning -----
rm(list = ls())
gc(reset = TRUE)

# renv -----
set.seed(42)

# paths/config -----
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)
source(here('config', 'global_config.R'))

# benchmarking -----
source(here('src', 'benchmark.R'))

# validation -----
source(here('shiny', 'common', 'validate.R'))
source(here('shiny', 'common', 'utils.R'))

# outlier detection (single source of truth) -----
source(here('shiny', 'wizard', 'outlier_operations.R'))

# project selection -----
# Allow overrides from environment for wizard
project <- Sys.getenv('PROJECT_NAME', unset = DEFAULT_PROJECT)

# Load config using centralized function (eliminates duplication)
source(here('shiny', 'wizard', 'config_operations.R'))
config_params <- loadConfigForPipeline(project)
list2env(config_params, envir = environment())

# functions -----
# Functions are loaded from shiny/common/utils.R

# Set up processing parameters
hours_eps <- as.numeric(Sys.getenv('HOURS_EPS', unset = HOURS_EPS))
use_logit <- as.logical(as.integer(Sys.getenv('LOGIT', unset = '0')))
logit_treat_inf <- as.logical(as.integer(Sys.getenv('LOGIT_TREAT_INF', unset = '1')))
use_iqr <- as.logical(as.integer(Sys.getenv('USE_IQR', unset = ifelse(isTRUE(USE_IQR), '1', '0'))))
tech_agg <- Sys.getenv('TECH_AGG', unset = TECH_AGG)
outlier_method <- Sys.getenv('OUTLIER_METHOD', unset = OUTLIER_METHOD)
outlier_clusters <- if (length(OUTLIER_CLUSTERS) > 0) OUTLIER_CLUSTERS else character(0)

# Start overall benchmark
benchmark_total <- startBenchmark("Processed Data Generation")

# Load raw data -----
benchmark_load <- startBenchmark("Raw Data Loading")
raw_file <- here('data', project, paste0(project, '_raw_merged_table.rds'))
if (!file.exists(raw_file)) {
  stop(sprintf("Raw data file not found: %s. Please run wizard/generate_raw.R first.", raw_file))
}
merged_table <- readRDS(raw_file)
endBenchmark(benchmark_load, sprintf("Loaded raw data: %d rows, %d columns", nrow(merged_table), ncol(merged_table)))

# percentage to logit transformation -----
merged_table <- applyLogitTransform(merged_table, isTRUE(use_logit), isTRUE(logit_treat_inf))

# outlier detection and treatment -----
benchmark_outliers <- startBenchmark("Outlier Detection")

# Get outlier variables from config or use defaults
outlier_variables <- if (exists('OUTLIER_DETECTION_VARIABLES') && length(OUTLIER_DETECTION_VARIABLES) > 0) {
  OUTLIER_DETECTION_VARIABLES
} else {
  # Default variables that commonly exist
  c('digital_biomass_mm3', 'height_mm', 'greenness_average')
}
# Only use variables that exist in the data
outlier_variables <- outlier_variables[outlier_variables %in% names(merged_table)]

# Get outlier factors from config or use defaults
outlier_factors <- if (exists('OUTLIER_DETECTION_FACTORS') && length(OUTLIER_DETECTION_FACTORS) > 0) {
  OUTLIER_DETECTION_FACTORS
} else {
  # Default factors that commonly exist
  c('dbscan_cluster', 'treatment')
}
outlier_factors <- outlier_factors[outlier_factors %in% names(merged_table)]

# Get detection method from config or use default
detection_method <- if (exists('OUTLIER_DETECTION_METHOD')) {
  OUTLIER_DETECTION_METHOD
} else {
  'iqr'  # Default to IQR method
}

if (length(outlier_variables) > 0) {
  # Use wizard's outlier detection as single source of truth
  if (detection_method == 'iqr' || detection_method == 'zscore') {
    # Apply per-cell outlier replacement using wizard functions
    merged_table <- replaceOutliersWithNA(
      data = merged_table,
      method = detection_method,
      variables = outlier_variables,
      factors = outlier_factors
    )
    
    cat(sprintf("Outlier detection completed using %s method on variables: %s\n", 
                detection_method, paste(outlier_variables, collapse = ', ')))
  } else {
    cat(sprintf("Outlier detection skipped: unsupported method '%s'\n", detection_method))
  }
} else {
  cat("Outlier detection skipped: no valid variables found\n")
}

endBenchmark(benchmark_outliers, "Outlier detection completed")

# aggregation for technical repetitions within time clusters -----
benchmark_agg <- startBenchmark("Technical Aggregation")

# Group by key factors and aggregate numeric columns
grouping_cols <- c('timestamp', 'unit', 'v_t_r', 'treatment', 'cultivar', 'dbscan_cluster', 'timestamp_group')
grouping_cols <- grouping_cols[grouping_cols %in% names(merged_table)]

cat(sprintf("Aggregating by columns: %s\n", paste(grouping_cols, collapse = ', ')))
cat(sprintf("Data before aggregation: %d rows, %d columns\n", nrow(merged_table), ncol(merged_table)))

agg_function <- if (tech_agg == 'median') median else mean

# Check if we actually need aggregation (multiple rows per group)
group_counts <- merged_table %>%
  dplyr::count(dplyr::across(dplyr::all_of(grouping_cols)))

if (any(group_counts$n > 1)) {
  cat("Multiple rows per group detected, performing aggregation...\n")
  merged_table <- merged_table %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), agg_function, na.rm = TRUE), .groups = 'drop')
} else {
  cat("No aggregation needed - single row per group\n")
}

endBenchmark(benchmark_agg, sprintf("Technical aggregation (%s): %d rows", tech_agg, nrow(merged_table)))

# medians for outliers within time clusters -----
benchmark_medians <- startBenchmark("Outlier Medians")

# Apply manual outlier cluster removal if specified
if (length(outlier_clusters) > 0) {
  outlier_table <- merged_table %>%
    dplyr::filter(dbscan_cluster %in% outlier_clusters)
  
  merged_table <- merged_table %>%
    dplyr::filter(!dbscan_cluster %in% outlier_clusters)
  
  cat(sprintf("Outlier clusters removed: %s (%d rows removed, %d remaining)\n", 
              paste(outlier_clusters, collapse = ', '), nrow(outlier_table), nrow(merged_table)))
} else {
  outlier_table <- data.frame()
}

endBenchmark(benchmark_medians, sprintf("Outlier handling: %d rows removed", nrow(outlier_table)))

# Data preprocessing: Fix single-level factors and ensure consistent data types -----
benchmark_preprocess <- startBenchmark("Data Preprocessing")
merged_table <- preprocessDataForAnalysis(merged_table)
endBenchmark(benchmark_preprocess, "Preprocessed data for analysis")

# exports -----
benchmark_export <- startBenchmark("Data Export")

final_table <- dplyr::bind_rows(
  merged_table %>% dplyr::mutate(outlier = FALSE),
  outlier_table %>% dplyr::mutate(outlier = TRUE)
)

dir.create(here('data', project), showWarnings = FALSE, recursive = TRUE)
saveRDS(final_table, file = here('data', project, paste0(project, '_merged_table.rds')))

# Save vector of groups (timestamp groups as POSIXct to prevent precision mismatch)
vector_of_groups <- if ('timestamp_group' %in% names(merged_table)) {
  sort(unique(merged_table$timestamp_group))  # Keep as POSIXct, not character
} else {
  # Fallback to character column names if timestamp_group not available
  merged_table %>%
    dplyr::select(where(is.character)) %>%
    colnames() %>%
    {
      .[!. %in% c("unit", "v_t_r", "cultivar")]
    }
}
saveRDS(vector_of_groups, file = here('data', project, paste0(project, '_vector_of_groups.rds')))

endBenchmark(benchmark_export, sprintf("Exported %d rows to .rds files", nrow(final_table)))

# End overall benchmark
endBenchmark(benchmark_total, sprintf("Complete pipeline for project %s", project))

# Performance check
cat("\n📊 Performance Summary:\n")
performance_check <- checkPerformanceTargets("Total Processing", 30)
cat(performance_check$message, "\n")

# Show benchmark summary
summary_stats <- getBenchmarkSummary()
if (!"message" %in% names(summary_stats)) {
  cat(sprintf("Total operations: %d, Average duration: %.3fs, Total time: %.3fs\n", 
              summary_stats$total_operations, summary_stats$avg_duration, summary_stats$total_duration))
}

remove(list = ls())
gc(reset = TRUE)

