# =============================================================================
# PROCESSING PIPELINE
# =============================================================================

#' Run the full data processing pipeline
#' @param merged_table Input data table
#' @param vector_of_groups Groups vector
#' @param config_params Configuration parameters
#' @return List with processed data and results
runProcessingPipeline <- function(merged_table, vector_of_groups, config_params) {
  logEvent('INFO', 'pipeline.start', list(
    nrow = nrow(merged_table),
    eps_hours = config_params$eps_hours,
    use_logit = isTRUE(config_params$use_logit),
    logit_treat_inf = isTRUE(config_params$logit_treat_inf),
    use_iqr = isTRUE(config_params$use_iqr),
    outlier_method = config_params$outlier_method,
    outlier_detection = config_params$outlier_detection,
    iqr_factors = paste(config_params$iqr_factors %||% character(0), collapse = ', '),
    tech_agg = config_params$tech_agg
  ))
  
  # Apply logit transform if enabled
  if (isTRUE(config_params$use_logit)) {
    merged_table <- applyLogitTransform(merged_table, TRUE, isTRUE(config_params$logit_treat_inf))
    logEvent('INFO', 'pipeline.logit.applied')
  }
  
  # Apply DBSCAN clustering (using centralized function)
  # Handle both field name variations
  eps <- config_params$eps_hours %||% config_params$hours_eps %||% 1
  merged_table <- addDBSCANClustersToData(merged_table, eps)
  if (is.null(merged_table)) {
    logError('pipeline.dbscan.failed')
    stop('DBSCAN clustering failed')
  }
  logEvent('INFO', 'pipeline.dbscan.applied', list(
    eps = eps,
    n_clusters = length(unique(merged_table$dbscan_cluster)),
    cluster_sizes = as.list(sort(table(merged_table$dbscan_cluster), decreasing = TRUE)[1:min(5, length(unique(merged_table$dbscan_cluster)))])
  ))
  
  # Apply outlier handling using unified strategy function
  if (isTRUE(config_params$use_iqr)) {
    factors <- config_params$iqr_factors %||% character(0)
    # Use all numeric columns for outlier detection (dynamic)
    variables <- names(merged_table)[sapply(merged_table, is.numeric)]
    
    logEvent('INFO', 'pipeline.outlier.start', list(
      strategy = config_params$outlier_method,
      method = config_params$outlier_detection,
      n_variables = length(variables),
      factors = paste(factors, collapse = ', ')
    ))
    
    # Apply unified outlier strategy
    result <- applyOutlierStrategy(
      data = merged_table,
      strategy = config_params$outlier_method,
      method = config_params$outlier_detection,
      variables = variables,
      factors = factors
    )
    
    merged_table <- result$data
    outlier_count <- result$outlier_count
    outlier_summary <- result$outlier_summary
    
    logEvent('INFO', 'pipeline.outlier.complete', list(
      strategy = config_params$outlier_method,
      outliers_detected = outlier_count,
      rows_after = nrow(merged_table)
    ))
  } else {
    merged_table$outlier <- FALSE
    outlier_count <- 0
    outlier_summary <- list()
  }
  
  # Apply technical aggregation
  if (nrow(merged_table) > 0) {
    # Check if required grouping columns exist
    required_cols <- c('unit', 'dbscan_cluster')
    missing_cols <- required_cols[!required_cols %in% names(merged_table)]
    if (length(missing_cols) > 0) {
      logError(paste('Missing required columns for aggregation:', paste(missing_cols, collapse = ', ')))
      stop(paste('Missing required columns for aggregation:', paste(missing_cols, collapse = ', ')))
    }
    
    agg_fun <- if (identical(tolower(config_params$tech_agg), 'mean')) mean else median
    before_rows <- nrow(merged_table)
    
    # Check if there are any numeric columns to aggregate
    numeric_cols <- names(merged_table)[sapply(merged_table, is.numeric)]
    if (length(numeric_cols) == 0) {
      logWarn('No numeric columns found for aggregation')
      # Just group and take first row of non-numeric data
      merged_table <- merged_table %>%
        dplyr::group_by(unit, dbscan_cluster) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    } else {
      merged_table <- merged_table %>%
        dplyr::group_by(unit, dbscan_cluster) %>%
        dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), function(x) agg_fun(x, na.rm = TRUE)), .groups = 'drop') %>%
        dplyr::left_join(
          merged_table %>%
            dplyr::select(-tidyselect::where(is.numeric)) %>%
            dplyr::group_by(unit, dbscan_cluster) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup(),
          by = c("unit", "dbscan_cluster")
        )
    }
    
    logEvent('INFO', 'pipeline.aggregation.applied', list(
      method = config_params$tech_agg,
      rows_before_agg = before_rows,
      rows_after_agg = nrow(merged_table)
    ))
  }
  
  # Remove DBSCAN outlier clusters (faulty timepoints)
  outlier_clusters <- config_params$outlier_clusters %||% character(0)
  if (length(outlier_clusters) > 0) {
    before_removal <- nrow(merged_table)
    outlier_clusters_num <- as.numeric(outlier_clusters)
    
    # Remove rows with selected outlier clusters
    merged_table <- merged_table %>%
      dplyr::filter(!as.numeric(as.character(dbscan_cluster)) %in% outlier_clusters_num)
    
    removed_rows <- before_removal - nrow(merged_table)
    
    logEvent('INFO', 'pipeline.dbscan_outliers.removed', list(
      outlier_clusters = paste(outlier_clusters, collapse = ', '),
      rows_removed = removed_rows,
      rows_remaining = nrow(merged_table)
    ))
  } else {
    logEvent('INFO', 'pipeline.dbscan_outliers.none_selected')
  }
  
  logEvent('INFO', 'pipeline.complete', list(final_rows = nrow(merged_table)))
  
  # Create updated vector_of_groups from processed data
  # This ensures timestamp groups match the actual data after DBSCAN cluster removal
  # Save as POSIXct to match the processed data exactly (prevents precision mismatch)
  updated_vector_of_groups <- if ('timestamp_group' %in% names(merged_table)) {
    sort(unique(merged_table$timestamp_group))  # Keep as POSIXct, not character
  } else {
    vector_of_groups  # Fallback to original if timestamp_group not available
  }
  
  return(list(
    merged_table = merged_table,
    vector_of_groups = updated_vector_of_groups,
    outlier_count = outlier_count,
    outlier_summary = outlier_summary
  ))
}

#' Save processed data files
#' @param processed_data List with merged_table and vector_of_groups
#' @param project_name Project name
#' @param export_name Export name suffix (optional)
#' @return List with file paths
saveProcessedData <- function(processed_data, project_name, export_name = '') {
  # Generate file names
  suffix <- if (nzchar(export_name)) paste0('_', export_name) else ''
  merged_file <- paste0(project_name, suffix, '_merged_table.rds')
  groups_file <- paste0(project_name, suffix, '_vector_of_groups.rds')
  
  # Save files in project folder (not shiny folder)
  project_dir <- here('data', project_name)
  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE)
  }
  
  merged_path <- here(project_dir, merged_file)
  groups_path <- here(project_dir, groups_file)
  
  saveRDS(processed_data$merged_table, merged_path)
  saveRDS(processed_data$vector_of_groups, groups_path)
  
  logEvent('INFO', 'pipeline.files.saved', list(
    merged_file = merged_file,
    groups_file = groups_file,
    merged_path = merged_path,
    groups_path = groups_path
  ))
  
  return(list(
    merged_file = merged_file,
    groups_file = groups_file,
    merged_path = merged_path,
    groups_path = groups_path
  ))
}
