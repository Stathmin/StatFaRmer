# StatFaRmer Minimal Deployment
# Pre-selected project version for ShinyApps.io

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
library(emmeans)
library(multcomp)
library(e1071)
library(lme4)
library(thematic)
library(jsonlite)

# Enable thematic for plot theming
thematic::thematic_shiny()

# Set options
options(shiny.reactlog = TRUE)
set.seed(42)

# Load here package for proper paths
library(here)

# Set the project option to pre-select project_NO3
options(statfarmer.project = "project_NO3")

# Load benchmarking and logging
source(here('src', 'benchmark.R'))
source(here('shiny', 'common', 'logger.R'))

# Load global config
source(here('config', 'global_config.R'))

# Source modular components
source(here('shiny', 'app', 'ui.R'))
source(here('shiny', 'app', 'server.R'))

# Run the application
shinyApp(ui = ui, server = server)
