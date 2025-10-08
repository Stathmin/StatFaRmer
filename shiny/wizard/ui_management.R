# =============================================================================
# WIZARD UI MANAGEMENT (SIMPLIFIED)
# =============================================================================
# Clean, simple UI update logic without excessive observers
#
# Replaces:
# - ui_controls.R: setupUIChoicesObserver() (324 lines → ~80 lines)
# - ui_controls.R: setupPendingObservers() (6 observers → 1 observer)
# - ui_controls.R: applyConfigSelections() (complex dual-format → centralized in config_operations.R)
# - ui_controls.R: setupValidationOutputs(), setupProjectDirectoryObserver(), setupValidationObserver()

#' Setup validation-related outputs
#' @param output Shiny output object
#' @param state State management object
setupValidationOutputs <- function(output, state) {
  output$validated <- shiny::renderText({ if (state$validated()) 'true' else 'false' })
  outputOptions(output, 'validated', suspendWhenHidden = FALSE)  # Needed for conditionalPanel
  output$validation_out <- renderText({ state$validation_log() })
}

#' Setup project directory population observer
#' @param session Shiny session
setupProjectDirectoryObserver <- function(session) {
  observe({
    dirs <- tryCatch({
      all <- list.dirs(here('data'), recursive = FALSE, full.names = FALSE)
      all[!grepl('^\\.', all)]
    }, error = function(e) character(0))
    updateSelectizeInput(session, 'input_dir', choices = sort(dirs), server = TRUE)
  })
}

#' Setup validation button observer
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session
#' @param state State management object
setupValidationObserver <- function(input, output, session, state) {
  observeEvent(input$validate_btn, {
    withBenchmark('Master:Validate', {
      wizardValidateProject(input$input_dir, state, session)
    })
  })
}

#' Setup UI choices population (SIMPLIFIED)
#' @param input Shiny input object
#' @param session Shiny session
#' @param state State management object
setupUIChoices <- function(input, session, state) {
  # Single reactive to track when UI needs update
  ui_trigger <- reactive({
    result <- list(
      project = input$input_dir,
      validated = state$validated(),
      use_logit = input$use_logit,
      logit_treat_inf = input$logit_treat_inf
    )
    logEvent('DEBUG', 'ui_trigger.reactive_calculated', list(
      project = result$project,
      validated = result$validated,
      use_logit = result$use_logit,
      logit_treat_inf = result$logit_treat_inf
    ))
    result
  })
  
  # Note: UI trigger observer removed to prevent reactive loops
  # All UI updates are now handled in validation_manager.R
}

#' Update DBSCAN clusters (simplified)
#' @param proj Project name
#' @param eps Epsilon parameter
#' @param use_logit Use logit transform
#' @param logit_treat_inf Treat infinite values
#' @param session Shiny session
#' @param state State management object
updateDBSCANClustersSimple <- function(proj, eps, use_logit, logit_treat_inf, session, state) {
  merged <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'dbscan.update')
  if (is.null(merged)) return()
  
  merged <- applyLogitTransform(merged, isTRUE(use_logit), isTRUE(logit_treat_inf))
  updateDBSCANClustersInUI(merged, eps %||% 1, session, state)
}

