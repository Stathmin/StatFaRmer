# StatFaRmer Master Wizard
# Entry point for data validation and preprocessing

library(shiny)
library(shinyjs)

# Check if we should launch the main app instead
project_file <- here::here('logs', 'launch_project.txt')
if (file.exists(project_file)) {
  # Read project name and clean up file
  project_name <- readLines(project_file, warn = FALSE)
  unlink(project_file)
  
  # Set the project option for the main app
  options(statfarmer.project = project_name)
  
  # Launch main app
  source(here::here('shiny', 'app.R'), local = TRUE)
} else {
  # Launch wizard
  source(here::here('shiny', 'wizard', 'master_ui.R'), local = TRUE)
  source(here::here('shiny', 'wizard', 'master_server.R'), local = TRUE)
  
  # Create wizard app
  shinyApp(ui = master_ui, server = master_server)
}

