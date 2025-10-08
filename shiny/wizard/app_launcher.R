# =============================================================================
# WIZARD APP LAUNCHER
# =============================================================================

#' Setup app launcher
#' 
#' @param input Shiny input object
#' @param session Shiny session object
#' @param state State management object
setupAppLauncher <- function(input, session, state) {
  
  observeEvent(input$launch_app_btn, {
    proj <- input$input_dir %||% ''
    if (!nzchar(proj)) return()
    
    # Simple deduplication: skip if already launching
    if (isTRUE(state$launching_app())) {
      logEvent('INFO', 'wizard.launch.skipped', list(reason = 'already_launching'))
      return()
    }
    
    state$launching_app(TRUE)
    removeModal()
    
    # Launch main app by reloading session
    logEvent('INFO', 'master.launch.app', list(project = proj))
    
    # Show notification to user
    showNotification(
      paste("Launching StatFaRmer for project:", proj),
      type = "message",
      duration = 2
    )
    
    # Set the project for the main app using a file
    project_file <- here::here('logs', 'launch_project.txt')
    writeLines(proj, project_file)
    
    # Reload the session to switch to main app
    session$reload()
  })
}
