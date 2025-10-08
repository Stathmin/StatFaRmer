# StatFaRmer Master Wizard (Iteration 4.1)

if (!requireNamespace('shiny', quietly = TRUE)) install.packages('shiny')
if (!requireNamespace('here', quietly = TRUE)) install.packages('here')
if (!requireNamespace('plotly', quietly = TRUE)) install.packages('plotly')
if (!requireNamespace('jsonlite', quietly = TRUE)) install.packages('jsonlite')
if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr')
if (!requireNamespace('stringi', quietly = TRUE)) install.packages('stringi')
library(shiny)
library(here)
library(plotly)
library(jsonlite)
library(dplyr)
library(stringi)

# Ensure logging/benchmarking available
source(here('src', 'benchmark.R'))
source(here('shiny', 'common', 'logger.R'))

# Load wizard UI/server
source(here('shiny', 'wizard', 'master_ui.R'))
source(here('shiny', 'wizard', 'master_server.R'))

masterApp <- function(host = '0.0.0.0', port = 3840) {
  shiny::shinyApp(ui = master_ui, server = master_server) |>
    shiny::runApp(host = host, port = port, launch.browser = FALSE)
}

if (identical(environment(), globalenv())) {
  # Allow direct `Rscript shiny/master_app.R`
  masterApp()
}




