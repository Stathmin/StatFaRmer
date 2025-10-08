# =============================================================================
# PLOT RENDERING (CONSOLIDATED)
# =============================================================================
# All plot functions in one place (merged from plot_helpers.R)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Prepare data for plotting
#' @param df Data frame
#' @param yvar Y variable name
#' @param input Shiny input object
#' @return Prepared data frame with display flags
prepareDataForPlot <- function(df, yvar, input) {
  if (!'outlier' %in% names(df)) {
    return(df)
  }
  
  # Combine automatic outliers and manually selected clusters
  selected_clusters <- input$outlier_clusters %||% character(0)
  if (length(selected_clusters) > 0) {
    sc_num <- suppressWarnings(as.integer(selected_clusters))
    df$.__manual_outlier__ <- as.integer(as.character(df$dbscan_cluster)) %in% sc_num
  } else {
    df$.__manual_outlier__ <- FALSE
  }
  
  # Robust vectorized flags (treat NA as FALSE)
  auto_flag <- ifelse(is.na(df$outlier), FALSE, as.logical(df$outlier))
  manual_flag <- ifelse(is.na(df$.__manual_outlier__), FALSE, as.logical(df$.__manual_outlier__))
  df$.__display_outlier__ <- auto_flag | manual_flag
  
  # Create legend labels
  df$.__legend_label__ <- ifelse(df$.__display_outlier__, 'Outlier', paste0('Cluster ', as.character(df$dbscan_cluster)))
  df$.__shape_label__ <- ifelse(df$.__display_outlier__, 'Outlier', 'Normal')
  
  return(df)
}

#' Create color scheme for plot
#' @param clusters Unique cluster values
#' @return Named vector of colors
createPlotColors <- function(clusters) {
  n_clusters <- length(clusters)
  cluster_colors <- rep(c('#000000', '#808080'), length.out = n_clusters)  # Black, gray alternating
  names(cluster_colors) <- as.character(clusters)
  
  # Add outlier color and rename clusters for legend
  legend_colors <- c('Outlier' = '#d62728')
  cluster_names <- paste0('Cluster ', names(cluster_colors))
  names(cluster_colors) <- cluster_names
  legend_colors <- c(legend_colors, cluster_colors)
  
  return(legend_colors)
}

#' Create base ggplot with outlier highlighting
#' @param df Data frame
#' @param yvar Y variable name
#' @param legend_colors Named color vector
#' @return ggplot object
createBasePlot <- function(df, yvar, legend_colors) {
  if ('.__legend_label__' %in% names(df)) {
    # Plot with outlier highlighting
    p <- ggplot2::ggplot(df, ggplot2::aes(x = timestamp, y = .data[[yvar]])) +
      ggplot2::geom_point(ggplot2::aes(color = .__legend_label__, shape = .__shape_label__), alpha = 0.7, size = 1.8) +
      ggplot2::scale_color_manual(values = legend_colors, breaks = names(legend_colors), name = 'Point Type') +
      ggplot2::scale_shape_manual(values = c('Normal' = 16, 'Outlier' = 17), name = 'Outlier') +
      ggplot2::labs(x = 'timestamp', y = yvar) +
      getStatfarmerTheme()
  } else {
    # Fallback plot without outlier highlighting
    clusters <- sort(unique(df$dbscan_cluster))
    cluster_colors <- rep(c('#000000', '#808080'), length.out = length(clusters))
    names(cluster_colors) <- as.character(clusters)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = timestamp, y = .data[[yvar]])) +
      ggplot2::geom_point(ggplot2::aes(color = as.factor(dbscan_cluster)), alpha = 0.6, size = 1) +
      ggplot2::scale_color_manual(values = cluster_colors) +
      ggplot2::labs(x = 'timestamp', y = yvar, color = 'DBSCAN Cluster') +
      getStatfarmerTheme()
  }
  
  return(p)
}

#' Apply outlier styling to plotly object
#' @param plotly_p Plotly object
#' @param input Shiny input object
#' @return Styled plotly object
applyOutlierStyling <- function(plotly_p, input) {
  outlier_clusters <- input$outlier_clusters %||% character(0)
  if (length(outlier_clusters) == 0) {
    return(plotly_p)
  }
  
  # Find traces that correspond to selected clusters
  cluster_traces <- which(sapply(plotly_p$x$data, function(x) {
    trace_name <- x$name
    if (is.null(trace_name)) return(FALSE)
    cluster_num <- gsub("Cluster ", "", trace_name)
    return(cluster_num %in% outlier_clusters)
  }))
  
  # Highlight selected clusters with red color
  if (length(cluster_traces) > 0) {
    for (trace_idx in cluster_traces) {
      plotly_p <- plotly_p %>%
        plotly::style(
          marker = list(size = 3, color = '#d62728', symbol = 'diamond'),
          traces = trace_idx
        )
    }
  }
  
  return(plotly_p)
}

#' Create empty plot with message
#' @param message Message to display
#' @return Plotly object
createEmptyPlot <- function(message = "No data available. Please select a project and variable.") {
  empty_plot <- ggplot2::ggplot() + 
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = message, size = 6) +
    getStatfarmerTheme() +
    ggplot2::theme_void()
  return(plotly::ggplotly(empty_plot))
}

# =============================================================================
# MAIN RENDERING FUNCTION
# =============================================================================

#' Render DBSCAN plot with outlier highlighting
#' @param plot_data List with data and yvar
#' @param input Shiny input object
#' @return Plotly object
renderDBSCANPlot <- function(plot_data, input) {
  # Handle null or empty data
  if (is.null(plot_data)) {
    return(createEmptyPlot("No data available. Please select a project and variable."))
  }
  
  df <- plot_data$data
  yvar <- plot_data$yvar
  
  if (nrow(df) == 0) {
    return(createEmptyPlot("No data points after filtering."))
  }
  
  # Prepare data with outlier flags and labels
  df <- prepareDataForPlot(df, yvar, input)
  
  # Create color scheme
  clusters <- sort(unique(df$dbscan_cluster))
  legend_colors <- createPlotColors(clusters)
  
  # Create base plot
  p <- createBasePlot(df, yvar, legend_colors)
  
  # Add facet if formula is valid
  ff <- tryCatch(as.formula(input$facet_formula), error = function(e) NULL)
  if (!is.null(ff)) {
    p <- p + ggplot2::facet_grid(ff)
  }
  
  # Convert to plotly
  plotly_p <- plotly::ggplotly(p, tooltip = c('x', 'y', 'colour', 'shape'))
  
  # Apply outlier styling
  plotly_p <- applyOutlierStyling(plotly_p, input)
  
  plotly_p
}
