# StatFaRmer Server Logic
# Iteration 3: Modular Shiny

# Load required libraries
library(shiny)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(broom)
library(flextable)
library(moments)
library(e1071)
library(viridis)
library(svglite)
library(DT)
library(glue)
library(lme4)
library(cowplot)
library(thematic)

# Load here package for proper paths
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)

# Enable thematic for plot theming
thematic::thematic_shiny()

# Load utility functions
source(here('shiny', 'common', 'utils.R'))
source(here('shiny', 'common', 'stats.R'))
source(here('shiny', 'common', 'plotting.R'))

# Ensure benchmarking and logging are available in server/env
source(here('src', 'benchmark.R'))
source(here('shiny', 'common', 'logger.R'))

# Load reactive modules
source(here('shiny', 'app', 'reactive_data.R'))
source(here('shiny', 'app', 'reactive_ui.R'))
source(here('shiny', 'app', 'reactive_outputs.R'))
source(here('shiny', 'app', 'reactive_downloads.R'))

# Load scientific modules
source(here('shiny', 'app', 'module_effect_sizes.R'))
source(here('shiny', 'app', 'module_growth_summaries.R'))

# Get available projects
available_projects <- getAvailableProjects()

# Default project honoring environment override and statfarmer.project option
default_project <- if (!is.null(getOption('statfarmer.project'))) {
  getOption('statfarmer.project')
} else {
  getDefaultProject(available_projects)
}

# Load default project data
data_list <- tryCatch({
  loadProjectData(default_project)
}, error = function(e) {
  # Fallback to old method if project-specific files don't exist
  loadShinyData()
})

merged_table <- data_list$merged_table
vector_of_groups <- data_list$vector_of_groups

# =============================================================================
# SERVER FUNCTION
# =============================================================================

server <- function(input, output, session) {
  set.seed(42)
  
  # Reactive values for combined inputs
  combined_inputs <- reactiveValues()
  
  # =============================================================================
  # REACTIVE DATA MANAGEMENT
  # =============================================================================
  
  # Create reactive data functions
  reactive_data <- createReactiveData(input, session, combined_inputs)
  projectData <- reactive_data$projectData
  projectUISelections <- reactive_data$projectUISelections
  filteredData <- reactive_data$filteredData

  # Project config JSON preview
  output$project_config_json <- renderText({
    sel <- input$selected_project %||% NA_character_
    if (is.na(sel) || !nzchar(sel)) return('No project selected')
    cfg_path <- here('data', sel, 'config.json')
    if (!file.exists(cfg_path)) return('config.json not found in project folder')
    tryCatch({
      paste(readLines(cfg_path, warn = FALSE), collapse = '\n')
    }, error = function(e) paste('Failed to read config.json:', e$message))
  })
  
  # =============================================================================
  # REACTIVE UI MANAGEMENT
  # =============================================================================
  
  # Create reactive UI functions
  reactive_ui <- createReactiveUI(input, output, session, combined_inputs, projectUISelections, projectData)
  
  # Set up UI observers
  reactive_ui$updateProjectSelections()
  reactive_ui$updateFactorLevels(projectData)
  reactive_ui$updateCultivars(projectData)
  reactive_ui$updateCombinedInputsInit()
  reactive_ui$updateCombinedInputsSubmit()
  
  # =============================================================================
  # REACTIVE OUTPUTS
  # =============================================================================
  
  # Create reactive output functions
  reactive_outputs <- createReactiveOutputs(input, output, session, combined_inputs, filteredData, projectData)
  analysisResults <- reactive_outputs$analysisResults
  
  # =============================================================================
  # REACTIVE DOWNLOADS
  # =============================================================================
  
  # Create reactive download functions
  createReactiveDownloads(input, output, session, combined_inputs, filteredData, analysisResults)
  
  # =============================================================================
  # SCIENTIFIC MODULES
  # =============================================================================
  
  # Effect Sizes Module
  effectSizesServer("effect_sizes", filteredData, analysisResults, combined_inputs)
  
  # Growth Summaries Module
  growthSummariesServer("growth_summaries", filteredData, combined_inputs, projectData)
}