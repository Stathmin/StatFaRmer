# StatFaRmer Shiny Application - Modular Version
# Iteration 3: Modular Shiny

# Load required libraries
library(shiny)
library(tidyverse)
library(stringr)
library(ggplot2)
library(multcompView)
library(broom)
library(flextable)
library(moments)
library(bslib)
library(viridis)
library(svglite)
library(shinyWidgets)
library(DT)
library(glue)
library(cowplot)

# Set options
options(shiny.reactlog = TRUE)
set.seed(42)

# Load here package for proper paths
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)

# Load benchmarking and logging
source(here('src', 'benchmark.R'))
source(here('shiny', 'common', 'logger.R'))

# Load global config
source(here('config', 'global_config.R'))

# Source modular components
source(here('shiny', 'app', 'ui.R'))
source(here('shiny', 'app', 'server.R'))

# Run the application
shiny::runApp(shinyApp(ui = ui, server = server),
              host = SHINY_HOST,
              port = SHINY_PORT,
              launch.browser = FALSE)
