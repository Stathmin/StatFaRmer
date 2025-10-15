# Raw Data Generation for Wizard
# This generates the raw data files that the wizard uses for consistent logit behavior

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
source(here('shiny', 'common', 'logger.R'))
source(here('shiny', 'common', 'validate.R'))
source(here('shiny', 'common', 'utils.R'))

# Function to generate raw data for a specific project
generateRawData <- function(project_name) {
  project <- project_name

# Load config using centralized function (eliminates duplication)
source(here('shiny', 'wizard', 'config_operations.R'))
config_params <- loadConfigForPipeline(project)
list2env(config_params, envir = environment())

# functions -----
# Functions are loaded from shiny/common/utils.R

# Project data validation
cat("🔍 Validating project data:", project, "\n")
validation_result <- validateProject(here('data', project))

if (!validation_result$valid) {
  cat("❌ Validation failed!\n")
  cat(formatValidationResult(validation_result), "\n")
  
  # Log validation errors
  logValidationResult(validation_result, project)
  
  # Graceful degradation - continue with warning
  cat("⚠️ Continuing execution with warning...\n")
} else {
  cat("✅ Validation passed successfully\n")
  logValidationResult(validation_result, project)
}

# Set up processing parameters
hours_eps <- as.numeric(Sys.getenv('HOURS_EPS', unset = HOURS_EPS))
use_logit <- as.logical(as.integer(Sys.getenv('LOGIT', unset = '0')))
logit_treat_inf <- as.logical(as.integer(Sys.getenv('LOGIT_TREAT_INF', unset = '1')))

# Set up file paths
project_dir <- here('data', project)
phenospex_file <- Sys.glob(paste0(project_dir, '/*_data.zip'))
unit_file_1 <- Sys.glob(paste0(project_dir, '/*_handmade.csv'))
unit_file_2 <- Sys.glob(paste0(project_dir, '/*_translation.csv'))
groups_file <- Sys.glob(paste0(project_dir, '/groups.xlsx'))

# Check if files exist
if (length(phenospex_file) == 0) {
  stop("No data zip file found in project directory")
}
if (length(unit_file_1) == 0) {
  stop("No handmade CSV file found in project directory")
}
if (length(unit_file_2) == 0) {
  stop("No translation CSV file found in project directory")
}
if (length(groups_file) == 0) {
  stop("No groups Excel file found in project directory")
}

# Start overall benchmark
benchmark_total <- startBenchmark("Raw Data Generation")

# import planteye table -----
benchmark_data_load <- startBenchmark("Data Loading")
planteye_table <- readr::read_csv(phenospex_file) %>%
  janitor::clean_names(.)

remove(phenospex_file)
endBenchmark(benchmark_data_load, sprintf("Loaded %d rows, %d columns", nrow(planteye_table), ncol(planteye_table)))

# aggregation with dbscan -----
benchmark_dbscan <- startBenchmark("DBSCAN Clustering")

dbscan_cluster <- tibble::as_tibble_col(
  forcats::as_factor(
    computeDbscanClustersFromTimestamps(planteye_table$timestamp, hours_eps)
  ), column_name = 'dbscan_cluster'
)

checkmate::assert_true(dbscan_cluster %>%
                     dplyr::n_distinct() > 1)

planteye_table <- planteye_table %>%
  dplyr::bind_cols(dbscan_cluster)

## legacy-compatible timestamp grouping for UI/facets
planteye_table <- planteye_table %>%
  dplyr::arrange(timestamp) %>%
  dplyr::group_by(dbscan_cluster) %>%
  dplyr::mutate(timestamp_group = mean(timestamp)) %>%
  dplyr::ungroup()

num_clusters <- length(unique(planteye_table$dbscan_cluster))
eps_value <- hours_eps
remove(list = c('dbscan_cluster', 'hours_eps'))
endBenchmark(benchmark_dbscan, sprintf("Created %d clusters with eps=%.1f hours", 
                                     num_clusters, eps_value))

# removal of rows with all observations at zero -----
planteye_table <- planteye_table %>%
  dplyr::filter(rowSums(dplyr::select(., where(is.numeric)) == 0, na.rm = TRUE) <
                  ncol(dplyr::select(., where(is.numeric))))

# import unit data -----
benchmark_unit <- startBenchmark("Unit Data Import")
unit_table_1 <- readr::read_csv(unit_file_1) %>%
  dplyr::mutate_all(as.character)
unit_table_2 <- readr::read_csv(unit_file_2) %>%
  dplyr::mutate_all(as.character)

unit_data <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = 'V.T.R') %>%
  janitor::clean_names(.)

remove(list = c('unit_file_1', 'unit_file_2', 'unit_table_1', 'unit_table_2'))
endBenchmark(benchmark_unit, sprintf("Loaded %d unit records", nrow(unit_data)))

# import groups table -----
benchmark_groups <- startBenchmark("Groups Import")
groups_table <- readxl::read_excel(groups_file) %>%
  janitor::clean_names(.)
remove(groups_file)
endBenchmark(benchmark_groups, sprintf("Loaded %d group records", nrow(groups_table)))

# table merge -----
benchmark_merge <- startBenchmark("Table Merge")

# Merge planteye with unit data
merged_table <- planteye_table %>%
  dplyr::left_join(unit_data, by = c('unit' = 't_x_y')) %>%
  dplyr::mutate(cultivar = as.character(cultivar)) %>%
  dplyr::left_join(groups_table %>% dplyr::mutate(cultivar = as.character(cultivar)), by = 'cultivar') %>%
  dplyr::select(-treatment.x) %>%
  dplyr::rename(treatment = treatment.y) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = 'UTC')) %>%
  dplyr::mutate(treatment = as.character(treatment))

# Clean up intermediate tables
remove(list = c('planteye_table', 'unit_data', 'groups_table'))

endBenchmark(benchmark_merge, sprintf("Merged tables: %d rows, %d columns", nrow(merged_table), ncol(merged_table)))

# Ensure required keys exist before aggregation -----
required_keys <- c('timestamp', 'unit', 'treatment', 'cultivar', 'dbscan_cluster')
missing_keys <- required_keys[!required_keys %in% names(merged_table)]

if (length(missing_keys) > 0) {
  cat("Available columns:", paste(names(merged_table), collapse = ', '), "\n")
  stop(sprintf("Missing required columns: %s", paste(missing_keys, collapse = ', ')))
}

# Data preprocessing: Fix single-level factors and ensure consistent data types -----
benchmark_preprocess <- startBenchmark("Data Preprocessing")
merged_table <- preprocessDataForAnalysis(merged_table)
endBenchmark(benchmark_preprocess, "Preprocessed data for analysis")

# table reordering -----
# Select available columns in order
available_cols <- c('timestamp', 'unit', 'treatment', 'cultivar', 'dbscan_cluster', 'timestamp_group')
available_cols <- available_cols[available_cols %in% names(merged_table)]
other_cols <- names(merged_table)[!names(merged_table) %in% available_cols]

merged_table <- merged_table %>%
  dplyr::select(all_of(available_cols), all_of(other_cols))

# Export raw data (before logit transform) -----
benchmark_export <- startBenchmark("Raw Data Export")

dir.create(here('data', project), showWarnings = FALSE, recursive = TRUE)

# Save raw merged table (with percent variables, no logit transform)
saveRDS(merged_table, file = here('data', project, paste0(project, '_raw_merged_table.rds')))

# Save raw vector of groups (timestamp groups as POSIXct to prevent precision mismatch)
raw_vector_of_groups <- if ('timestamp_group' %in% names(merged_table)) {
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
saveRDS(raw_vector_of_groups, file = here('data', project, paste0(project, '_raw_vector_of_groups.rds')))

endBenchmark(benchmark_export, sprintf("Exported raw data: %d rows to .rds files", nrow(merged_table)))

# End overall benchmark
endBenchmark(benchmark_total, sprintf("Raw data generation for project %s", project))

# Performance check
cat("\n📊 Raw Generation Performance Summary:\n")
performance_check <- checkPerformanceTargets("Raw Generation", 15)
cat(performance_check$message, "\n")

# Show benchmark summary
summary_stats <- getBenchmarkSummary()
if (!"message" %in% names(summary_stats)) {
  cat(sprintf("Total operations: %d, Average duration: %.3fs, Total time: %.3fs\n", 
              summary_stats$total_operations, summary_stats$avg_duration, summary_stats$total_duration))
}

  remove(list = ls())
  gc(reset = TRUE)
}

# Check if project_name variable is set and call the function
if (exists('project_name') && nzchar(project_name)) {
  cat("🚀 Starting raw data generation for project:", project_name, "\n")
  generateRawData(project_name)
  cat("✅ Raw data generation completed for project:", project_name, "\n")
} else {
  cat("ℹ️ project_name variable not set. Function defined but not executed.\n")
  cat("   To execute: project_name <- 'your_project' then source this file.\n")
}
