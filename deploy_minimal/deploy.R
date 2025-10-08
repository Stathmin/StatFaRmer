# StatFaRmer Minimal Deployment Script
# For ShinyApps.io deployment with pre-selected project

# Install required packages if not already installed
required_packages <- c(
  "shiny", "tidyverse", "ggplot2", "emmeans", "multcomp", 
  "multcompView", "broom", "flextable", "moments", "bslib", 
  "viridis", "svglite", "shinyWidgets", "DT", "glue", "cowplot", 
  "e1071", "lme4", "thematic", "jsonlite", "here"
)

# Install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load libraries
library(shiny)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(multcomp)
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
library(e1071)
library(lme4)
library(thematic)
library(jsonlite)
library(here)

# Enable thematic for plot theming
thematic::thematic_shiny()

# Set options
options(shiny.reactlog = TRUE)
set.seed(42)

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
