# StatFaRmer Global Configuration

# Paths
DATA_DIR <- here::here('data')
SHINY_DIR <- here::here('shiny')
LOGS_DIR <- here::here('logs')

# Cache & Debug
ENABLE_CACHE <- TRUE
ENABLE_DEBUG_CACHE <- FALSE

# Benchmark targets
TARGET_TOTAL_SECONDS <- 30

# Processing defaults
DEFAULT_PROJECT <- 'project_NO3'
HOURS_EPS <- 1
USE_IQR <- FALSE

# Shiny app defaults - SECURE by default
SHINY_HOST <- '127.0.0.1'  # Localhost only - secure default
SHINY_PORT <- 3839

# Validation
REQUIRED_FIELDS <- list(
  traitfinder = c('timestamp', 'unit'),
  metadata = c('V.T.R', 'Treatment', 'Cultivar'),
  translation = c('V.T.R', 'T:X:Y')
)


