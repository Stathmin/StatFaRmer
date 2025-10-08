library(shiny)
library(here)
library(plotly)
library(dbscan)
library(cowplot)
library(thematic)

# Enable thematic for plot theming
thematic::thematic_shiny()

# Source common modules
source(here('shiny', 'common', 'utils.R'))
source(here('shiny', 'common', 'validate.R'))
source(here('shiny', 'common', 'stats.R'))
source(here('shiny', 'common', 'logger.R'))
source(here('shiny', 'common', 'plotting.R'))
source(here('src', 'benchmark.R'))

# Source wizard operation modules
source(here('shiny', 'wizard', 'data_operations.R'))
source(here('shiny', 'wizard', 'clustering_operations.R'))
source(here('shiny', 'wizard', 'outlier_operations.R'))
source(here('shiny', 'wizard', 'config_operations.R'))

# Source remaining wizard modules
source(here('shiny', 'wizard', 'winsorizing.R'))
source(here('shiny', 'wizard', 'process_pipeline.R'))
source(here('shiny', 'wizard', 'plot_rendering.R'))
source(here('shiny', 'wizard', 'minimal_preprocess.R'))

# Source simplified wizard modules
source(here('shiny', 'wizard', 'state_management.R'))
source(here('shiny', 'wizard', 'ui_management.R'))
source(here('shiny', 'wizard', 'validation_manager.R'))
source(here('shiny', 'wizard', 'data_processor.R'))
source(here('shiny', 'wizard', 'pipeline_manager.R'))
source(here('shiny', 'wizard', 'app_launcher.R'))

# =============================================================================
# MASTER SERVER (REFACTORED)
# =============================================================================

master_server <- function(input, output, session) {
  logEvent('INFO', 'master.start', list())
  
  # Initialize state management
  state <- initWizardState()
  
  # Setup simplified UI management
  logEvent('INFO', 'master_server.setup_ui', list())
  setupValidationOutputs(output, state)
  setupProjectDirectoryObserver(session)
  setupValidationObserver(input, output, session, state)
  setupUIChoices(input, session, state)
  setupEpsilonObserver(input, session, state)
  setupDBSCANNumericObserver(input, session, state)
  setupLogitTransformObserver(input, session, state)
  # Add back UI observers but with loading flag to prevent loops
  setupUIStateObservers(input, session, state)
  
  # Setup pipeline execution
  setupPipelineExecution(input, output, state)
  
  # Setup app launcher
  setupAppLauncher(input, session, state)
  
  # Create reactive data and outputs
  plotData <- createPlotDataReactive(input, state)
  paramsReactive <- createParamsReactive(input, state)
  
  # Setup outputs
  output$dbscan_plot <- plotly::renderPlotly({
    renderDBSCANPlot(plotData(), input)
  })
  
  output$params_out <- renderText({
    paste(paramsReactive(), collapse = '\n')
  })
}
