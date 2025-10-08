# =============================================================================
# WIZARD VALIDATION MANAGEMENT
# =============================================================================

#' Handle project validation
#' 
#' Validates project folder, creates raw cache, and loads configuration.
#' Updates state variables based on validation results.
#' 
#' @param proj Project name
#' @param state State management object
#' @param session Shiny session
#' @return TRUE if validation successful, FALSE otherwise
wizardValidateProject <- function(proj, state, session) {
  if (!nzchar(proj)) {
    state$append_validation('Please select a project folder.')
    logWarn('master.validate: empty project selection')
    state$validated(FALSE)
    return(FALSE)
  }
  
  proj_path <- here('data', proj)
  if (!dir.exists(proj_path)) {
    logError(paste('Project folder does not exist:', proj_path))
    state$append_validation(paste('ERROR: folder not found -', proj_path))
    state$validated(FALSE)
    return(FALSE)
  }
  
  # Use real project validation
  state$validation_log("")
  state$append_validation(
    'Starting validation for ', proj, '\n',
    '- Project path: ', proj_path, '\n',
    '- Scanning files...'
  )
  
  res <- tryCatch(validateProject(proj_path), error = function(e) list(valid = FALSE, messages = e$message))
  
  if (is.list(res) && isTRUE(res$valid)) {
    state$append_validation('Validation passed.')
    
    # Summarize project files
    summarizeProjectFiles(proj, proj_path, state)
    
    # Create raw cache
    ok <- createRawCache(proj, state)
    
    # Set validation state (this will trigger UI to appear)
    state$validated(ok)
    
    # Load and apply configuration using delayed approach for UI timing
    # This must happen AFTER state$validated(ok) so conditionalPanel renders
    tryCatch({
      cfg <- loadProjectConfig(proj)
      if (!is.null(cfg)) {
        logEvent('INFO', 'validation.config_loaded', list(
          project = proj,
          has_pp = !is.null(cfg$processing_parameters),
          has_os = !is.null(cfg$outlier_settings)
        ))
        
        # Update state immediately (doesn't depend on UI)
        updateAppliedFromConfig(cfg, state)
        
        # Schedule UI updates after conditionalPanel has rendered
        if (requireNamespace('shinyjs', quietly = TRUE)) {
          logEvent('DEBUG', 'shinyjs.delay.starting', list(project = proj))
          shinyjs::delay(300, {
            logEvent('INFO', 'validation.applying_config_to_ui', list(project = proj))
                
                # Set loading flag to prevent observers from firing
                state$loading_config(TRUE)
                logEvent('DEBUG', 'loading_config.set', list(flag = TRUE))
                
                # Apply basic config if available
                if (!is.null(cfg)) {
                  applyBasicConfig(cfg, session)
                }
            
            # Always populate UI choices regardless of config availability
            merged <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'validation')
            if (!is.null(merged)) {
              # Apply logit transform if config is available, otherwise use defaults
              use_logit <- if (!is.null(cfg)) isTRUE(cfg$processing_parameters$use_logit) else FALSE
              logit_treat_inf <- if (!is.null(cfg)) isTRUE(cfg$processing_parameters$logit_treat_inf) else TRUE
              merged <- applyLogitTransform(merged, use_logit, logit_treat_inf)
              
              # Update UI choices (factors and numerics)
              factors <- names(merged)[sapply(merged, function(x) is.character(x) || is.factor(x))]
              numerics <- names(merged)[sapply(merged, is.numeric)]
              if ('dbscan_cluster' %in% names(merged) && !'dbscan_cluster' %in% factors) {
                factors <- c(factors, 'dbscan_cluster')
              }
              
              # Set iqr_factors choices and selection (empty if no config)
              selected_factors <- if (!is.null(cfg) && !is.null(cfg$processing_parameters$iqr_factors) && length(cfg$processing_parameters$iqr_factors) > 0) {
                as.character(cfg$processing_parameters$iqr_factors)
              } else {
                character(0)  # Empty selection when no config
              }
              updateSelectizeInput(session, 'iqr_factors', choices = factors, selected = selected_factors, server = TRUE)
              # Also set the state for iqr_factors
              if (!is.null(state) && !is.null(state$applied)) {
                state$applied$iqr_factors <- selected_factors
              }
              
              # Set other config values directly (use defaults if no config)
              if (!is.null(state) && !is.null(state$applied)) {
                state$applied$use_iqr <- if (!is.null(cfg)) isTRUE(cfg$processing_parameters$use_iqr) else FALSE
                state$applied$outlier_method <- if (!is.null(cfg)) cfg$outlier_settings$outlier_method %||% 'remove' else 'remove'
                state$applied$outlier_detection <- if (!is.null(cfg)) cfg$outlier_settings$outlier_detection %||% 'iqr' else 'iqr'
                state$applied$tech_agg <- if (!is.null(cfg)) tolower(cfg$processing_parameters$tech_agg %||% 'median') else 'median'
              }
              updateSelectizeInput(session, 'dbscan_numeric', choices = numerics, server = TRUE)
              
              # Auto-select dbscan_numeric with logit transform logic
              current_selection <- session$input$dbscan_numeric
              new_selection <- NULL
              
              if (!is.null(current_selection) && current_selection %in% numerics) {
                # Keep current selection if it's still available
                new_selection <- current_selection
              } else if (!is.null(current_selection)) {
                # Try to find logit equivalent
                if (grepl('_percent$', current_selection)) {
                  logit_equivalent <- gsub('_percent$', '_logit', current_selection)
                  if (logit_equivalent %in% numerics) {
                    new_selection <- logit_equivalent
                  }
                } else if (grepl('_logit$', current_selection)) {
                  percent_equivalent <- gsub('_logit$', '_percent', current_selection)
                  if (percent_equivalent %in% numerics) {
                    new_selection <- percent_equivalent
                  }
                }
              }
              
              # Default to first numeric column if no selection and no config
              if (is.null(new_selection)) {
                if (!is.null(cfg) && 'digital_biomass_mm3' %in% numerics) {
                  new_selection <- 'digital_biomass_mm3'
                } else if (length(numerics) > 0) {
                  new_selection <- numerics[1]  # First numeric column
                }
              }
              
              if (!is.null(new_selection)) {
                updateSelectizeInput(session, 'dbscan_numeric', selected = new_selection)
                # Also set the state for plot data reactive
                if (!is.null(state) && !is.null(state$applied)) {
                  state$applied$dbscan_numeric <- new_selection
                }
              }
              
              # KISS: Directly calculate clusters and set both choices and selection
              eps_hours <- if (!is.null(cfg)) cfg$processing_parameters$hours_eps %||% 1 else 1
              tmp <- addDBSCANClustersToData(merged, eps_hours)
              if (!is.null(tmp)) {
                new_clusters <- sort(unique(as.integer(as.character(tmp$dbscan_cluster))))
                cluster_choices <- as.character(new_clusters)
                names(cluster_choices) <- paste0('Cluster ', new_clusters)
                
                # Set both choices and selection in one call (empty if no config)
                selected_clusters <- if (!is.null(cfg) && !is.null(cfg$outlier_settings$outlier_clusters) && length(cfg$outlier_settings$outlier_clusters) > 0) {
                  # Handle both single and multiple outlier_clusters cases
                  # Both string and list cases are handled the same way with as.character()
                  as.character(cfg$outlier_settings$outlier_clusters)
                } else {
                  character(0)  # Empty selection when no config
                }
                
                updateSelectizeInput(session, 'outlier_clusters',
                                  choices = cluster_choices,
                                  selected = selected_clusters,
                                  server = TRUE)
                
                logEvent('INFO', 'config.outlier_clusters.direct', list(
                  choices_count = length(cluster_choices),
                  selected = paste(selected_clusters, collapse = ', '),
                  has_config = !is.null(cfg),
                  loading_config = isTRUE(state$loading_config()),
                  eps_hours = eps_hours
                ))
              }
            }
            
                # Apply other config selections only if config exists
                if (!is.null(cfg)) {
                  applyConfigSelections(proj, session, state)
                }
                
                # Clear loading flag to allow observers to fire again
                state$loading_config(FALSE)
                logEvent('DEBUG', 'loading_config.cleared', list(flag = FALSE))
              })
        } else {
          logWarn('shinyjs not available for delayed config application')
        }
      } else {
        logEvent('INFO', 'validation.no_config', list(project = proj))
        
        # Still populate UI choices even without config
        if (requireNamespace('shinyjs', quietly = TRUE)) {
          logEvent('DEBUG', 'shinyjs.delay.starting_no_config', list(project = proj))
          shinyjs::delay(300, {
            logEvent('INFO', 'validation.populating_ui_no_config', list(project = proj))
            
            # Set loading flag to prevent observers from firing
            state$loading_config(TRUE)
            logEvent('DEBUG', 'loading_config.set', list(flag = TRUE))
            
            # Load data and populate UI choices with defaults
            merged <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'validation_no_config')
            if (!is.null(merged)) {
              # Use default logit transform settings
              merged <- applyLogitTransform(merged, FALSE, TRUE)
              
              # Update UI choices (factors and numerics)
              factors <- names(merged)[sapply(merged, function(x) is.character(x) || is.factor(x))]
              numerics <- names(merged)[sapply(merged, is.numeric)]
              if ('dbscan_cluster' %in% names(merged) && !'dbscan_cluster' %in% factors) {
                factors <- c(factors, 'dbscan_cluster')
              }
              
              # Set iqr_factors choices with empty selection
              updateSelectizeInput(session, 'iqr_factors', choices = factors, selected = character(0), server = TRUE)
              if (!is.null(state) && !is.null(state$applied)) {
                state$applied$iqr_factors <- character(0)
              }
              
              # Set other default values
              if (!is.null(state) && !is.null(state$applied)) {
                state$applied$use_iqr <- FALSE
                state$applied$outlier_method <- 'remove'
                state$applied$outlier_detection <- 'iqr'
                state$applied$tech_agg <- 'median'
              }
              
              # Set dbscan_numeric choices and select first numeric column
              updateSelectizeInput(session, 'dbscan_numeric', choices = numerics, server = TRUE)
              if (length(numerics) > 0) {
                first_numeric <- numerics[1]
                updateSelectizeInput(session, 'dbscan_numeric', selected = first_numeric)
                if (!is.null(state) && !is.null(state$applied)) {
                  state$applied$dbscan_numeric <- first_numeric
                }
              }
              
              # Calculate DBSCAN clusters with default epsilon
              tmp <- addDBSCANClustersToData(merged, 1)  # Default epsilon = 1
              if (!is.null(tmp)) {
                new_clusters <- sort(unique(as.integer(as.character(tmp$dbscan_cluster))))
                cluster_choices <- as.character(new_clusters)
                names(cluster_choices) <- paste0('Cluster ', new_clusters)
                
                # Set choices with empty selection
                updateSelectizeInput(session, 'outlier_clusters',
                                  choices = cluster_choices,
                                  selected = character(0),
                                  server = TRUE)
                
                logEvent('INFO', 'config.outlier_clusters.no_config', list(
                  choices_count = length(cluster_choices),
                  selected = 'none'
                ))
              }
            }
            
            # Clear loading flag to allow observers to fire again
            state$loading_config(FALSE)
            logEvent('DEBUG', 'loading_config.cleared', list(flag = FALSE))
          })
        }
      }
    }, error = function(e) {
      logError(paste('validation.config_error:', e$message))
    })
    
    logEvent('INFO', 'master.validated', list(project = proj))
    return(ok)
  } else {
    msg <- if (is.list(res) && !is.null(res$messages)) paste(res$messages, collapse = '\n') else 'Validation failed.'
    state$append_validation(msg)
    state$validated(FALSE)
    return(FALSE)
  }
}

#' Summarize project files and cache status
#' 
#' @param proj Project name
#' @param proj_path Project path
#' @param state State management object
summarizeProjectFiles <- function(proj, proj_path, state) {
  # Find project files
  zip_files <- list.files(proj_path, pattern = "*_data\\.zip$", full.names = TRUE)
  handmade_files <- list.files(proj_path, pattern = "*_handmade\\.csv$", full.names = TRUE)
  translation_files <- list.files(proj_path, pattern = "*_translation\\.csv$", full.names = TRUE)
  groups_files <- list.files(proj_path, pattern = "groups\\.xlsx$", full.names = TRUE)
  
  # Cache paths
  cache_merged <- here('data', proj, paste0(proj, '_merged_table.rds'))
  cache_groups <- here('data', proj, paste0(proj, '_vector_of_groups.rds'))
  raw_merged <- here('data', proj, paste0(proj, '_raw_merged_table.rds'))
  raw_groups <- here('data', proj, paste0(proj, '_raw_vector_of_groups.rds'))
  
  # Preview cache sizes if available
  merged_rows <- NA_integer_; merged_cols <- NA_integer_
  if (file.exists(cache_merged)) {
    prev <- tryCatch(readRDS(cache_merged), error = function(e) NULL)
    if (!is.null(prev)) { 
      merged_rows <- nrow(prev); merged_cols <- ncol(prev) 
    }
  }
  
  # Log summary in batch
  state$append_validation(
    'Preparing raw cache...\n',
    '- Files found: ', length(zip_files), ' zip, ', length(handmade_files), ' handmade, ', 
    length(translation_files), ' translation, ', length(groups_files), ' groups\n',
    '- Cache present: merged=', file.exists(cache_merged), ', groups=', file.exists(cache_groups), '\n',
    if (!is.na(merged_rows)) paste0('- Cache sizes: rows=', merged_rows, ', cols=', merged_cols, '\n') else '',
    '- Raw cache present: merged=', file.exists(raw_merged), ', groups=', file.exists(raw_groups)
  )
}

#' Create raw cache for project
#' 
#' @param proj Project name
#' @param state State management object
#' @return TRUE if successful, FALSE otherwise
createRawCache <- function(proj, state) {
  state$generating_raw(TRUE)
  
  ok <- TRUE
  tryCatch({
    raw_paths <- minimalPreprocess(proj)
    logEvent('INFO', 'master.raw.available', list(project = proj, raw = raw_paths))
    state$append_validation(
      sprintf('Raw cache created for %s (rows=%s, cols=%s)\n', proj, as.character(raw_paths$rows), as.character(raw_paths$cols)),
      '- Paths:\n',
      sprintf('  %s\n', raw_paths$raw_merged),
      sprintf('  %s', raw_paths$raw_groups)
    )
  }, error = function(e) {
    ok <<- FALSE
    logError(paste('master.raw.failed', e$message))
    state$append_validation('Validation passed, but raw cache failed: ', e$message)
  })
  
  state$generating_raw(FALSE)
  return(ok)
}

#' Setup epsilon change observer for DBSCAN recalculation
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupEpsilonObserver <- function(input, session, state) {
  observeEvent(input$eps_hours, {
    logEvent('DEBUG', 'epsilon_observer.fired', list(
      eps_hours = input$eps_hours,
      validated = isTRUE(state$validated()),
      loading_config = isTRUE(state$loading_config())
    ))
    if (!isTRUE(state$validated())) return()
    if (isTRUE(state$loading_config())) {
      logEvent('DEBUG', 'epsilon_observer.skipped_loading_config')
      return()  # Skip during config loading
    }
    proj <- input$input_dir %||% ''
    if (!nzchar(proj)) return()
    
          # Clear DBSCAN cache when epsilon changes
          state$dbscan_cache$data <- NULL
          state$dbscan_cache$eps <- NULL
          state$dbscan_cache$clusters <- NULL
          
          # Recalculate DBSCAN clusters when epsilon changes
          merged <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'epsilon_change')
          if (!is.null(merged)) {
            merged <- applyLogitTransform(merged, isTRUE(input$use_logit), isTRUE(input$logit_treat_inf))
            tmp <- addDBSCANClustersToData(merged, input$eps_hours %||% 1)
            if (!is.null(tmp)) {
              # Recalculate timestamp groups from new DBSCAN clusters
              # Each DBSCAN cluster gets one timestamp group (mean timestamp)
              tmp <- tmp %>%
                dplyr::arrange(timestamp) %>%
                dplyr::group_by(dbscan_cluster) %>%
                dplyr::mutate(timestamp_group = mean(timestamp)) %>%
                dplyr::ungroup()
              
              new_clusters <- sort(unique(as.integer(as.character(tmp$dbscan_cluster))))
              cluster_choices <- as.character(new_clusters)
              names(cluster_choices) <- paste0('Cluster ', new_clusters)
              
              # Preserve current selection if it still exists
              current_selected <- input$outlier_clusters %||% character(0)
              valid_selected <- current_selected[current_selected %in% cluster_choices]
              
              updateSelectizeInput(session, 'outlier_clusters', 
                                choices = cluster_choices, 
                                selected = valid_selected,
                                server = TRUE)
              
              # Update the DBSCAN cache with the new data including recalculated timestamp groups
              state$dbscan_cache$data <- tmp
              state$dbscan_cache$eps <- input$eps_hours %||% 1
              state$dbscan_cache$clusters <- new_clusters
              
              logEvent('INFO', 'epsilon_observer.timestamp_groups_recalculated', list(
                eps = input$eps_hours %||% 1,
                n_clusters = length(new_clusters),
                n_timestamp_groups = length(unique(tmp$timestamp_group))
              ))
            }
          }
  }, ignoreInit = TRUE)
}

#' Setup dbscan_numeric observer to update state
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupDBSCANNumericObserver <- function(input, session, state) {
  observeEvent(input$dbscan_numeric, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$dbscan_numeric)) {
      state$applied$dbscan_numeric <- input$dbscan_numeric
    }
  }, ignoreInit = TRUE)
}

#' Setup logit transform observer to clear DBSCAN cache
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupLogitTransformObserver <- function(input, session, state) {
  observeEvent(c(input$use_logit, input$logit_treat_inf), {
    # Clear DBSCAN cache when logit transform settings change
    state$dbscan_cache$data <- NULL
    state$dbscan_cache$eps <- NULL
    state$dbscan_cache$clusters <- NULL
  }, ignoreInit = TRUE)
}

#' Setup iqr_factors observer to update state
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupIqrFactorsObserver <- function(input, session, state) {
  observeEvent(input$iqr_factors, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$iqr_factors)) {
      state$applied$iqr_factors <- input$iqr_factors
    }
  }, ignoreInit = TRUE)
}

#' Setup observers for all UI controls that need to update state
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupAllUIObservers <- function(input, session, state) {
  # Tech aggregation observer
  observeEvent(input$tech_agg, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$tech_agg)) {
      state$applied$tech_agg <- input$tech_agg
    }
  }, ignoreInit = TRUE)
  
  # Outlier detection method observer
  observeEvent(input$outlier_detection, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$outlier_detection)) {
      state$applied$outlier_detection <- input$outlier_detection
    }
  }, ignoreInit = TRUE)
  
  # Outlier method observer
  observeEvent(input$outlier_method, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$outlier_method)) {
      state$applied$outlier_method <- input$outlier_method
    }
  }, ignoreInit = TRUE)
  
  # Use IQR observer
  observeEvent(input$use_iqr, {
    if (!is.null(state) && !is.null(state$applied) && !is.null(input$use_iqr)) {
      state$applied$use_iqr <- input$use_iqr
    }
  }, ignoreInit = TRUE)
  
  # Outlier clusters observer
  observeEvent(input$outlier_clusters, {
    # This is handled by the epsilon observer, but we need to ensure it's saved to config
    # The config saving happens in pipeline_manager.R
  }, ignoreInit = TRUE)
}

#' Setup UI state observers with loading flag to prevent loops
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupUIStateObservers <- function(input, session, state) {
  # IQR factors observer
  observeEvent(input$iqr_factors, {
    logEvent('DEBUG', 'iqr_factors.observer.fired', list(
      loading_config = state$loading_config(),
      input_value = paste(input$iqr_factors, collapse = ', '),
      state_exists = !is.null(state) && !is.null(state$applied)
    ))
    if (!isTRUE(state$loading_config()) && !is.null(state) && !is.null(state$applied)) {
      state$applied$iqr_factors <- input$iqr_factors
      logEvent('DEBUG', 'iqr_factors.state.updated', list(
        new_value = paste(input$iqr_factors, collapse = ', ')
      ))
    }
  }, ignoreInit = TRUE)
  
  # Tech aggregation observer
  observeEvent(input$tech_agg, {
    if (!isTRUE(state$loading_config()) && !is.null(state) && !is.null(state$applied)) {
      state$applied$tech_agg <- input$tech_agg
    }
  }, ignoreInit = TRUE)
  
  # Outlier detection method observer
  observeEvent(input$outlier_detection, {
    if (!isTRUE(state$loading_config()) && !is.null(state) && !is.null(state$applied)) {
      state$applied$outlier_detection <- input$outlier_detection
    }
  }, ignoreInit = TRUE)
  
  # Outlier method observer
  observeEvent(input$outlier_method, {
    if (!isTRUE(state$loading_config()) && !is.null(state) && !is.null(state$applied)) {
      state$applied$outlier_method <- input$outlier_method
    }
  }, ignoreInit = TRUE)
  
  # Use IQR observer
  observeEvent(input$use_iqr, {
    if (!isTRUE(state$loading_config()) && !is.null(state) && !is.null(state$applied)) {
      state$applied$use_iqr <- input$use_iqr
    }
  }, ignoreInit = TRUE)
}

