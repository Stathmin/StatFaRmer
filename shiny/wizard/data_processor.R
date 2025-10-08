# =============================================================================
# WIZARD DATA PROCESSING
# =============================================================================

#' Create reactive plot data
#' 
#' Creates reactive data for plotting with DBSCAN clusters and outlier detection.
#' Handles data loading, transformations, and filtering.
#' 
#' @param input Shiny input object
#' @param state State management object
#' @return Reactive function that returns processed data
createPlotDataReactive <- function(input, state) {
  reactive({
    if (isTRUE(state$generating_raw())) return(NULL)
    
    proj <- input$input_dir %||% ''
    yvar <- if (!is.null(state) && !is.null(state$applied)) state$applied$dbscan_numeric %||% '' else ''
    eps <- input$eps_hours %||% 1
    
    if (!nzchar(proj) || !nzchar(yvar)) return(NULL)
    
    # Load project data
    df <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'plotData')
    if (is.null(df)) return(NULL)
    
    # Apply logit transform if enabled (optimized for plotting - only transform selected variable)
    if (isTRUE(input$use_logit)) {
      df <- applyLogitTransformForPlot(df, yvar, isTRUE(input$logit_treat_inf))
    }
    
    # Add DBSCAN clusters (with caching)
    df <- getCachedDBSCANData(df, eps, state)
    if (is.null(df)) return(NULL)
    
    # Check if yvar exists after transformations, or use logit version if available
    plot_yvar <- yvar
    if (isTRUE(input$use_logit) && grepl('_percent$', yvar)) {
      logit_yvar <- stringi::stri_replace_all_fixed(yvar, pattern = '_percent', replacement = '_logit')
      if (logit_yvar %in% names(df)) {
        plot_yvar <- logit_yvar
      }
    }
    
    if (!plot_yvar %in% names(df)) {
      logError(paste('plotData: Y variable not found:', plot_yvar))
      return(NULL)
    }
    
    # Compute outlier flags (use original yvar for outlier detection)
    detection_method <- input$outlier_detection %||% 'iqr'
    grouping_factors <- if (!is.null(state) && !is.null(state$applied)) state$applied$iqr_factors %||% character(0) else character(0)
    if (length(grouping_factors) > 0) {
      grouping_factors <- grouping_factors[grouping_factors %in% names(df)]
    }
    df <- computeOutlierFlags(df, yvar, detection_method, grouping_factors, add_outlier_col = TRUE)
    
    # Apply outlier filtering if requested
    df <- applyOutlierFiltering(df, input)
    
    return(list(data = df, yvar = plot_yvar))
  })
}
#' Apply outlier filtering to data
#' 
#' @param df Data frame
#' @param input Shiny input object
#' @return Filtered data frame
applyOutlierFiltering <- function(df, input) {
  if (isTRUE(input$hide_outliers_plot)) {
    # Hide automatic outliers (from outlier detection)
    if ('outlier' %in% names(df)) {
      df <- df[!df$outlier, ]
    }
    
    # Also hide manually selected clusters
    outlier_clusters <- input$outlier_clusters %||% character(0)
    if (length(outlier_clusters) > 0) {
      # Convert to numeric for comparison
      outlier_clusters_num <- as.numeric(outlier_clusters)
      df <- df[!as.numeric(as.character(df$dbscan_cluster)) %in% outlier_clusters_num, ]
    }
  }
  
  return(df)
}

#' Create parameters display reactive
#' 
#' @param input Shiny input object
#' @param state State management object
#' @return Reactive function that returns parameter display text
createParamsReactive <- function(input, state) {
  reactive({
    proj <- input$input_dir %||% ''
    outlier_clusters <- input$outlier_clusters %||% character(0)
    iqr_factors <- if (!is.null(state) && !is.null(state$applied)) state$applied$iqr_factors %||% character(0) else character(0)
    hidden_points <- if (length(outlier_clusters) > 0) {
      paste('Hidden clusters:', paste(outlier_clusters, collapse = ', '))
    } else {
      'No hidden clusters'
    }
    
    list(
      paste('=== Processing Configuration ==='),
      paste('project =', proj),
      paste('eps_hours =', input$eps_hours),
      paste('use_logit =', input$use_logit),
      paste('logit_treat_inf =', input$logit_treat_inf),
      paste('use_iqr =', if (!is.null(state) && !is.null(state$applied)) state$applied$use_iqr else FALSE),
      paste('outlier_method =', if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_method else 'remove'),
      paste('outlier_detection =', if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr'),
      paste('iqr_factors =', paste(iqr_factors, collapse = ', ')),
      paste('tech_agg =', if (!is.null(state) && !is.null(state$applied)) state$applied$tech_agg else 'median'),
      paste('facet_formula =', input$facet_formula),
      paste(''),
      paste('=== Outlier Settings ==='),
      paste('outlier_clusters =', paste(outlier_clusters, collapse = ', ')),
      paste(hidden_points),
      paste('hide_outliers_plot =', input$hide_outliers_plot),
      paste(''),
      paste('=== Export ==='),
      paste('export_name =', input$export_name)
    )
  })
}

#' Apply logit transform for plotting (optimized - only selected variable)
#' @param data Data frame
#' @param yvar Variable to transform
#' @param treat_inf Handle infinities
#' @return Data frame with logit-transformed variable
applyLogitTransformForPlot <- function(data, yvar, treat_inf = TRUE) {
  # Only transform the selected variable if it's a percentage column
  if (!yvar %in% names(data)) return(data)
  
  # Check if this is a percentage column that needs transformation
  if (!grepl('_percent$', yvar)) return(data)
  
  # Create logit version of the variable name
  logit_var <- stringi::stri_replace_all_fixed(yvar, pattern = '_percent', replacement = '_logit')
  
  # Apply transformation only to the selected variable
  data[[logit_var]] <- stats::qlogis(dplyr::case_when(
    (data[[yvar]] >= 0) & (data[[yvar]] <= 1.00) ~ data[[yvar]],
    (data[[yvar]] < 0) & (data[[yvar]] >= -0.01) ~ 0,
    (data[[yvar]] > 1) & (data[[yvar]] <= 1.01) ~ 1,
    .default = NA
  ))
  
  # Handle infinities if requested
  if (isTRUE(treat_inf)) {
    finite_vals <- data[[logit_var]][is.finite(data[[logit_var]])]
    if (length(finite_vals) > 0) {
      min_finite <- min(finite_vals, na.rm = TRUE)
      max_finite <- max(finite_vals, na.rm = TRUE)
      data[[logit_var]] <- pmax(pmin(data[[logit_var]], max_finite), min_finite)
    }
  }
  
  return(data)
}

#' Get cached DBSCAN data to avoid recalculating clusters
#' @param data Data frame
#' @param eps Epsilon parameter for DBSCAN
#' @param state State management object
#' @return Data frame with DBSCAN clusters
getCachedDBSCANData <- function(data, eps, state) {
  # Check if we can use cached result
  if (!is.null(state$dbscan_cache$data) && 
      !is.null(state$dbscan_cache$eps) && 
      !is.null(state$dbscan_cache$clusters) &&
      state$dbscan_cache$eps == eps &&
      identical(state$dbscan_cache$data, data)) {
    # Return cached result
    return(state$dbscan_cache$clusters)
  }
  
  # Calculate new DBSCAN clusters
  result <- addDBSCANClustersToData(data, eps)
  if (!is.null(result)) {
    # Cache the result
    state$dbscan_cache$data <- data
    state$dbscan_cache$eps <- eps
    state$dbscan_cache$clusters <- result
  }
  
  return(result)
}

