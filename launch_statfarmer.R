#!/usr/bin/env Rscript
# StatFaRmer Launcher
# Launches Master Wizard, then optionally launches main app with validated project

suppressPackageStartupMessages({
  library(shiny)
  library(here)
})

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  ____  _        _   _____      ____  __  __            \n")
cat(" / ___|| |_ __ _| |_|  ___|_ _ |  _ \\|  \\/  | ___ _ __ \n")
cat(" \\___ \\| __/ _` | __| |_ / _` || |_) | |\\/| |/ _ \\ '__|\n")
cat("  ___) | || (_| | |_|  _| (_| ||  _ <| |  | |  __/ |   \n")
cat(" |____/ \\__\\__,_|\\__|_|  \\__,_||_| \\_\\_|  |_|\\___|_|   \n")
cat("                                                         \n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")
cat("Statistical Analysis Platform for Plant Phenotyping\n")
cat("\n")
cat("WORKFLOW:\n")
cat("  1. Master Wizard: Validate & process project data\n")
cat("  2. Main App: Statistical analysis & visualization\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")

# Function to launch wizard
launch_wizard <- function() {
  cat("Launching Master Wizard...\n\n")
  tryCatch({
    result <- runApp(here('shiny', 'wizard'), launch.browser = TRUE)
    
    # Check if wizard returned a transition request
    if (is.list(result) && !is.null(result$action) && result$action == "launch_app") {
      cat("\n═══════════════════════════════════════════════════════════════\n")
      cat("Wizard completed successfully! Transitioning to main app...\n")
      cat("═══════════════════════════════════════════════════════════════\n\n")
      
      # Launch main app with the project from wizard
      launch_app(result$project)
    }
  }, error = function(e) {
    cat("ERROR launching wizard:", e$message, "\n")
    cat("Make sure shiny/wizard/ directory exists\n")
  })
}

# Function to launch main app with project
launch_app <- function(project = NULL) {
  if (is.null(project)) {
    cat("Launching main app (no project pre-selected)...\n\n")
    runApp(here('shiny', 'app'), launch.browser = TRUE)
  } else {
    cat(paste0("Launching main app with project: ", project, "\n\n"))
    # Pass project via environment for app default selection
    Sys.setenv(SELECTED_PROJECT = project)
    runApp(here('shiny', 'app'), launch.browser = TRUE)
  }
}

# Check if project was provided as command line arg
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  action <- args[1]
  
  if (action == "wizard") {
    launch_wizard()
  } else if (action == "app") {
    project <- if (length(args) > 1) args[2] else NULL
    launch_app(project)
  } else {
    cat("Unknown action:", action, "\n")
    cat("Usage:\n")
    cat("  Rscript launch_statfarmer.R wizard          # Launch Master Wizard\n")
    cat("  Rscript launch_statfarmer.R app [project]   # Launch main app\n")
  }
} else {
  # Interactive menu
  cat("What would you like to do?\n\n")
  cat("  [1] Launch Master Wizard (data validation & processing)\n")
  cat("  [2] Launch Main App (statistical analysis)\n")
  cat("  [Q] Quit\n\n")
  
  choice <- readline(prompt = "Enter your choice: ")
  
  if (tolower(choice) == "1" || tolower(choice) == "wizard") {
    launch_wizard()
  } else if (tolower(choice) == "2" || tolower(choice) == "app") {
    # Check for available projects
    data_dir <- here('data')
    if (dir.exists(data_dir)) {
      projects <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
      projects <- projects[startsWith(projects, 'project_')]
      
      if (length(projects) > 0) {
        cat("\nAvailable projects:\n")
        for (i in seq_along(projects)) {
          # Check if RDS exists
          rds_path <- here('data', projects[i], 'merged_table.rds')
          status <- if (file.exists(rds_path)) "✓" else "✗"
          cat(sprintf("  [%d] %s %s\n", i, status, projects[i]))
        }
        cat("\nEnter project number (or press Enter to select in app): ")
        proj_choice <- readline()
        
        if (nzchar(proj_choice) && !is.na(as.integer(proj_choice))) {
          idx <- as.integer(proj_choice)
          if (idx >= 1 && idx <= length(projects)) {
            launch_app(projects[idx])
          } else {
            launch_app()
          }
        } else {
          launch_app()
        }
      } else {
        cat("\nNo projects found. Launch wizard first to create one.\n")
        cat("Launching main app anyway...\n")
        launch_app()
      }
    } else {
      launch_app()
    }
  } else {
    cat("Goodbye!\n")
  }
}

