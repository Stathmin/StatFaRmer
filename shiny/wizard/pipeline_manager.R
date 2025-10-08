# =============================================================================
# WIZARD PIPELINE MANAGEMENT
# =============================================================================

#' Setup pipeline execution
#' 
#' Sets up the run button observer and handles the complete processing pipeline.
#' Manages outlier detection, config generation, and data processing.
#' 
#' @param input Shiny input object
#' @param output Shiny output object
#' @param state State management object
setupPipelineExecution <- function(input, output, state) {
  
  # Setup run button observer
  observeEvent(input$run_btn, {
    withBenchmark('Master:Run', {
      proj <- input$input_dir %||% ''
      
      # Validate project selection
      if (!nzchar(proj)) {
        output$log_tail <- renderText('Please select a project folder before running.')
        logWarn('master.run: empty project selection')
        return()
      }
      
      # Run outlier detection if enabled
      outlier_results <- runOutlierDetectionIfEnabled(proj, state, input)
      
      # Generate project configuration
      if (!generateProjectConfiguration(proj, input, state, outlier_results, output)) return()
      
      # Run processing pipeline
      if (!executeProcessingPipeline(proj, input, output, state)) return()
      
      # Show completion modal
      showCompletionModal(proj)
      
      # Update log output
      tail <- tryCatch(readLines(here('logs','app.log')), error = function(e) character(0))
      output$log_tail <- renderText(paste(tail[max(1, length(tail)-50):length(tail)], collapse = '\n'))
    })
  })
}

#' Run outlier detection if enabled
#' 
#' @param proj Project name
#' @param state State management object
#' @param input Shiny input object
#' @return Outlier results or NULL
runOutlierDetectionIfEnabled <- function(proj, state, input) {
  outlier_results <- NULL
  
  if (!is.null(state) && !is.null(state$applied) && isTRUE(state$applied$use_iqr)) {
    # Load project data for outlier detection
    avail <- tryCatch(loadProjectData(proj), error = function(e) NULL)
    if (!is.null(avail)) {
      merged_table <- avail$merged_table
      factors <- if (!is.null(state) && !is.null(state$applied)) state$applied$iqr_factors %||% character(0) else character(0)
      
      # If dbscan_cluster is requested as a grouping factor, compute it here
      if ('dbscan_cluster' %in% factors && 'timestamp' %in% names(merged_table) && requireNamespace('dbscan', quietly = TRUE)) {
        eps <- input$eps_hours %||% 1
        merged_table <- addDBSCANClustersToData(merged_table, eps)
        if (is.null(merged_table)) {
          logError('pipeline_manager.dbscan.failed')
        }
      }
      
      # Use all numeric columns for outlier detection (dynamic)
      variables <- names(merged_table)[sapply(merged_table, is.numeric)]
      
      # For config generation, we only need a lightweight check - not full outlier processing
      # The actual outlier handling will be done in the pipeline
      outlier_results <- list(
        outliers = data.frame(),  # Empty - not needed for config
        summary = list(),         # Empty - not needed for config  
        method = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr',
        factors = factors,
        variables = variables,
        total_outliers = 0        # Placeholder - actual count will be computed in pipeline
      )
      logEvent('INFO', 'master.outlier.detection', list(
        method = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr',
        factors = paste(factors, collapse = ', '),
        total_outliers = "computed_in_pipeline"
      ))
    }
  }
  
  return(outlier_results)
}

#' Generate project configuration
#' 
#' @param proj Project name
#' @param input Shiny input object
#' @param state State management object
#' @param outlier_results Outlier detection results
#' @param output Shiny output object
#' @return TRUE if successful, FALSE otherwise
generateProjectConfiguration <- function(proj, input, state, outlier_results, output) {
  config_params <- list(
    eps_hours = input$eps_hours,
    use_iqr = if (!is.null(state) && !is.null(state$applied)) state$applied$use_iqr else FALSE,
    use_logit = input$use_logit,
    logit_treat_inf = input$logit_treat_inf,
    tech_agg = if (!is.null(state) && !is.null(state$applied)) state$applied$tech_agg else 'median',
    iqr_factors = if (!is.null(state) && !is.null(state$applied)) state$applied$iqr_factors %||% character(0) else character(0),
    outlier_clusters = input$outlier_clusters %||% character(0),
    outlier_method = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_method else 'remove',
    outlier_detection = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr',
    facet_formula = input$facet_formula,
    outlier_results = outlier_results
  )
  
  config_success <- generateProjectConfig(proj, config_params)
  if (!config_success) {
    output$log_tail <- renderText('Failed to generate project configuration.')
    logError('master.run: config generation failed')
    return(FALSE)
  }
  
  return(TRUE)
}

#' Execute processing pipeline
#' 
#' @param proj Project name
#' @param input Shiny input object
#' @param output Shiny output object
#' @return TRUE if successful, FALSE otherwise
executeProcessingPipeline <- function(proj, input, output, state) {
  # Load project data (use raw data for processing)
  merged_table <- tryCatch(loadProjectDataWithFallback(proj, prefer = 'raw', context = 'pipeline'), error = function(e) {
    logError(paste('master.run.load_failed', e$message))
    stop(e)
  })
  
  if (is.null(merged_table)) {
    output$log_tail <- renderText('Failed to load project data.')
    return(FALSE)
  }
  
  # Load groups vector from raw data
  raw_groups_path <- here('data', proj, paste0(proj, '_raw_vector_of_groups.rds'))
  vector_of_groups <- tryCatch(readRDS(raw_groups_path), error = function(e) {
    logError(paste('master.run.groups_failed', e$message))
    stop(e)
  })
  
  # Run processing pipeline: writes treated RDS under data/<project>/
  tryCatch({
    # Get config params from the config that was already generated
    config_params <- list(
      eps_hours = input$eps_hours,
      use_iqr = if (!is.null(state) && !is.null(state$applied)) state$applied$use_iqr else FALSE,
      use_logit = input$use_logit,
      logit_treat_inf = input$logit_treat_inf,
      tech_agg = if (!is.null(state) && !is.null(state$applied)) state$applied$tech_agg else 'median',
      iqr_factors = if (!is.null(state) && !is.null(state$applied)) state$applied$iqr_factors %||% character(0) else character(0),
      outlier_clusters = input$outlier_clusters %||% character(0),
      outlier_method = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_method else 'remove',
      outlier_detection = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr',
      facet_formula = input$facet_formula
    )
    
    processed_data <- runProcessingPipeline(merged_table, vector_of_groups, config_params)
    files_saved <- saveProcessedData(processed_data, proj, input$export_name)
    
    # Update config with actual outlier count from pipeline
    if (!is.null(processed_data$outlier_count) && processed_data$outlier_count > 0) {
      updateConfigWithOutlierCount(proj, processed_data$outlier_count, processed_data$outlier_summary)
    }
    
    logEvent('INFO', 'master.run.complete', list(
      project = proj,
      files = paste(names(files_saved), collapse = ', '),
      outliers_detected = processed_data$outlier_count %||% 0
    ))
    
    return(TRUE)
    
  }, error = function(e) {
    logError(paste('Pipeline processing failed:', e$message))
    output$log_tail <- renderText(paste('ERROR during processing:', e$message))
    return(FALSE)
  })
}

#' Show completion modal
#' 
#' @param proj Project name
showCompletionModal <- function(proj) {
  showModal(modalDialog(
    title = 'Processing complete',
    easyClose = TRUE,
    footer = tagList(
      modalButton('Close'),
      actionButton('launch_app_btn', 'Launch Main App')
    ),
    div(
      p('Data processed successfully for project: ', strong(proj)),
      p('You can launch the main app now. It will open with this project pre-selected.'),
      p('Note: a short delay is applied to prevent redundant launches.')
    )
  ))
}
