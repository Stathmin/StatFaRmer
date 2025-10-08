# =============================================================================
# WIZARD CONFIGURATION OPERATIONS
# =============================================================================
# Configuration loading, saving, and management functions

#' Load project configuration
#' @param proj Project name
#' @return Configuration list or NULL if not found
loadProjectConfig <- function(proj) {
  config_path <- here('data', proj, 'config.json')
  
  if (!file.exists(config_path)) {
    logEvent('INFO', 'config.loading.not_found', list(project = proj))
    return(NULL)
  }
  
  tryCatch({
    config_data <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
    logEvent('INFO', 'config.loading.found', list(
      project = proj,
      path = config_path,
      has_processing_params = !is.null(config_data$processing_parameters),
      has_outlier_settings = !is.null(config_data$outlier_settings),
      names = paste(names(config_data), collapse = ", ")
    ))
    return(config_data)
  }, error = function(e) {
    logError(paste('config.loading.failed:', e$message))
    return(NULL)
  })
}

#' Load config for pipeline scripts (generate_raw.R, generate_processed.R)
#' Returns environment variables with defaults from config.json or global defaults
#' @param project Project name
#' @return Named list with all pipeline parameters
loadConfigForPipeline <- function(project) {
  cfg <- loadProjectConfig(project)
  if (is.null(cfg)) cfg <- list()  # Use empty list for defaults
  
  # Extract values with defaults
  pp <- cfg$processing_parameters %||% list()
  os <- cfg$outlier_settings %||% list()
  od <- cfg$outlier_detection %||% list()
  
  return(list(
    HOURS_EPS = pp$hours_eps %||% 1,
    USE_IQR = pp$use_iqr %||% FALSE,
    USE_LOGIT = pp$use_logit %||% FALSE,
    LOGIT_TREAT_INF = pp$logit_treat_inf %||% TRUE,
    TECH_AGG = pp$tech_agg %||% 'mean',
    IQR_FACTORS = pp$iqr_factors %||% character(0),
    FACET_FORMULA = pp$facet_formula %||% 'treatment ~ .',
    OUTLIER_CLUSTERS = os$outlier_clusters %||% character(0),
    OUTLIER_METHOD = os$outlier_method %||% 'iqr',
    OUTLIER_DETECTION_METHOD = od$outlier_detection_method %||% 'iqr',
    OUTLIER_DETECTION_FACTORS = od$outlier_detection_factors %||% character(0),
    OUTLIER_DETECTION_VARIABLES = od$outlier_detection_variables %||% character(0),
    TOTAL_OUTLIERS_DETECTED = od$total_outliers_detected %||% 0,
    OUTLIERS_BY_VARIABLE = od$outliers_by_variable %||% list()
  ))
}

#' Apply basic configuration values to UI
#' @param cfg Configuration object
#' @param session Shiny session
applyBasicConfig <- function(cfg, session) {
  pp <- cfg$processing_parameters
  os <- cfg$outlier_settings
  
  if (!is.null(pp$hours_eps)) updateNumericInput(session, 'eps_hours', value = as.numeric(pp$hours_eps))
  if (!is.null(pp$use_logit)) updateCheckboxInput(session, 'use_logit', value = isTRUE(pp$use_logit))
  if (!is.null(pp$logit_treat_inf)) updateCheckboxInput(session, 'logit_treat_inf', value = isTRUE(pp$logit_treat_inf))
  if (!is.null(pp$use_iqr)) updateCheckboxInput(session, 'use_iqr', value = isTRUE(pp$use_iqr))
  if (!is.null(pp$tech_agg)) updateRadioButtons(session, 'tech_agg', selected = tolower(pp$tech_agg))
  if (!is.null(os$outlier_method)) updateRadioButtons(session, 'outlier_method', selected = os$outlier_method)
  if (!is.null(os$outlier_detection)) updateRadioButtons(session, 'outlier_detection', selected = os$outlier_detection)
}

#' Update applied state from configuration
#' @param cfg Configuration object
#' @param state State management object
updateAppliedFromConfig <- function(cfg, state) {
  pp <- cfg$processing_parameters
  os <- cfg$outlier_settings
  
  if (!is.null(state) && !is.null(state$applied)) {
    state$applied$use_iqr <- isTRUE(pp$use_iqr)
    state$applied$outlier_method <- os$outlier_method %||% 'remove'
    state$applied$outlier_detection <- os$outlier_detection %||% 'iqr'
    state$applied$iqr_factors <- pp$iqr_factors %||% character(0)
    state$applied$tech_agg <- tolower(pp$tech_agg %||% 'median')
  }
  
  logEvent('INFO', 'config.applied.values', list(
    use_iqr = if (!is.null(state) && !is.null(state$applied)) state$applied$use_iqr else FALSE,
    outlier_method = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_method else 'remove',
    outlier_detection = if (!is.null(state) && !is.null(state$applied)) state$applied$outlier_detection else 'iqr',
    iqr_factors = if (!is.null(state) && !is.null(state$applied)) paste(state$applied$iqr_factors, collapse = ", ") else "",
    tech_agg = if (!is.null(state) && !is.null(state$applied)) state$applied$tech_agg else 'median'
  ))
}


#' Apply configuration selections to UI
#' @param proj Project name
#' @param session Shiny session
#' @param state State management object
applyConfigSelections <- function(proj, session, state) {
  cfg <- loadProjectConfig(proj)
  if (is.null(cfg)) {
    logEvent('WARN', 'config.selections.not_found', list(project = proj))
    return()
  }
  
  pp <- cfg$processing_parameters
  os <- cfg$outlier_settings
  
  # Note: iqr_factors is now handled directly in validation_manager.R after choices are populated
  
  # Apply facet formula
  if (!is.null(pp$facet_formula)) {
    updateTextInput(session, 'facet_formula', value = pp$facet_formula)
    logEvent('INFO', 'config.facet_formula.applied', list(selected = pp$facet_formula))
  }
  
  # Note: outlier_clusters is now handled directly in validation_manager.R after DBSCAN update
}

#' Generate project configuration file as JSON
#' @param project_name Project name
#' @param config_params List of configuration parameters
#' @return TRUE if successful, FALSE otherwise
generateProjectConfig <- function(project_name, config_params) {
  tryCatch({
    project_dir <- here('data', project_name)
    if (!dir.exists(project_dir)) {
      logError(paste('Project directory does not exist:', project_dir))
      return(FALSE)
    }

    config_path <- here(project_dir, 'config.json')

    # Create config data structure
    config_data <- list(
      metadata = list(
        generated_by = "StatFaRmer Master Wizard",
        generated_at = as.character(Sys.time()),
        project = project_name
      ),
      processing_parameters = list(
        hours_eps = config_params$eps_hours,
        use_iqr = isTRUE(config_params$use_iqr),
        use_logit = isTRUE(config_params$use_logit),
        logit_treat_inf = isTRUE(config_params$logit_treat_inf),
        tech_agg = config_params$tech_agg,
        iqr_factors = config_params$iqr_factors %||% character(0),
        facet_formula = config_params$facet_formula
      ),
      outlier_settings = list(
        outlier_clusters = config_params$outlier_clusters %||% character(0),
        outlier_method = config_params$outlier_method,
        outlier_detection = config_params$outlier_detection
      )
    )
    
    # Add outlier detection results if available
    if (!is.null(config_params$outlier_results)) {
      config_data$outlier_detection <- generateOutlierConfig(config_params$outlier_results)
    }

    # Write JSON config file
    jsonlite::write_json(config_data, config_path, pretty = TRUE, auto_unbox = TRUE)
    logEvent('INFO', 'master.config.generated', list(project = project_name, path = config_path))
    return(TRUE)
  }, error = function(e) {
    logError(paste('Failed to generate project config:', e$message))
    return(FALSE)
  })
}

#' Update config file with actual outlier count from pipeline
#' @param project_name Project name
#' @param outlier_count Total number of outliers detected
#' @param outlier_summary Summary of outliers by variable
#' @return TRUE if successful, FALSE otherwise
updateConfigWithOutlierCount <- function(project_name, outlier_count, outlier_summary) {
  tryCatch({
    config_path <- here('data', project_name, 'config.json')
    
    if (!file.exists(config_path)) {
      logWarn('Config file not found for outlier count update')
      return(FALSE)
    }
    
    # Read existing config
    config_data <- jsonlite::fromJSON(config_path)
    
    # Update outlier detection section
    if (is.null(config_data$outlier_detection)) {
      config_data$outlier_detection <- list()
    }
    
    config_data$outlier_detection$total_outliers_detected <- outlier_count
    config_data$outlier_detection$outliers_by_variable <- outlier_summary
    
    # Write updated config
    jsonlite::write_json(config_data, config_path, pretty = TRUE, auto_unbox = TRUE)
    logEvent('INFO', 'config.outlier_count.updated', list(
      project = project_name,
      total_outliers = outlier_count
    ))
    return(TRUE)
  }, error = function(e) {
    logError(paste('Failed to update config with outlier count:', e$message))
    return(FALSE)
  })
}
