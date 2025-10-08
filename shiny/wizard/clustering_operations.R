# =============================================================================
# WIZARD CLUSTERING OPERATIONS
# =============================================================================
# DBSCAN clustering and related operations
#
# Functions consolidated from:
# - data_processor.R: addDBSCANClusters()
# - ui_controls.R: updateDBSCANClusters(), applyStagedClusters()
# - common/utils.R: addDBSCANClustersToData() (moved here for better organization)

#' Add DBSCAN clusters to data (consolidated function)
#' @param data Data frame with timestamp column
#' @param eps Epsilon parameter for DBSCAN
#' @return Data frame with dbscan_cluster column added
addDBSCANClustersToData <- function(data, eps) {
  if (!requireNamespace('dbscan', quietly = TRUE)) {
    logError('DBSCAN: dbscan package not available')
    return(NULL)
  }
  
  if (!'timestamp' %in% names(data)) {
    logError('DBSCAN: timestamp column not found')
    return(NULL)
  }
  
  data$dbscan_cluster <- as.factor(computeDbscanClustersFromTimestamps(data$timestamp, eps))
  data <- renumberClustersByTime(data)
  
  # Recalculate timestamp groups from DBSCAN clusters
  # Each DBSCAN cluster gets one timestamp group (mean timestamp)
  data <- data %>%
    dplyr::arrange(timestamp) %>%
    dplyr::group_by(dbscan_cluster) %>%
    dplyr::mutate(timestamp_group = mean(timestamp)) %>%
    dplyr::ungroup()
  
  n_clusters <- length(unique(data$dbscan_cluster))
  n_timestamp_groups <- length(unique(data$timestamp_group))
  logEvent('INFO', 'dbscan.complete', list(
    eps = eps, 
    n_clusters = n_clusters,
    n_timestamp_groups = n_timestamp_groups
  ))
  
  return(data)
}

#' Update DBSCAN clusters in UI
#' @param merged Merged data table
#' @param eps Epsilon parameter
#' @param session Shiny session
updateDBSCANClustersInUI <- function(merged, eps, session, state = NULL) {
  # Use centralized DBSCAN function
  tmp <- addDBSCANClustersToData(merged, eps)
  if (is.null(tmp)) return()
  
  new_clusters <- sort(unique(as.integer(as.character(tmp$dbscan_cluster))))
  cluster_choices <- as.character(new_clusters)
  names(cluster_choices) <- paste0('Cluster ', new_clusters)
  
  logEvent('DEBUG', 'dbscan.ui_update', list(
    n_clusters = length(new_clusters),
    cluster_choices = paste(cluster_choices, collapse = ', '),
    eps = eps
  ))
  
  # Preserve currently selected values that still exist in new choices
  current_selected <- session$input$outlier_clusters %||% character(0)
  valid_selected <- current_selected[current_selected %in% cluster_choices]
  
  logEvent('DEBUG', 'dbscan.selectize_update', list(
    current_selected = paste(current_selected, collapse = ', '),
    valid_selected = paste(valid_selected, collapse = ', '),
    choices_count = length(cluster_choices)
  ))
  
  updateSelectizeInput(session, 'outlier_clusters', 
                      choices = cluster_choices, 
                      selected = valid_selected,
                      server = TRUE)
}
