# =============================================================================
# WIZARD STATE MANAGEMENT (SIMPLIFIED)
# =============================================================================

#' Initialize wizard state management
#' @return List containing essential reactive state variables
initWizardState <- function() {
  
  # Core state
  validated <- reactiveVal(FALSE)
  generating_raw <- reactiveVal(FALSE)
  validation_log <- reactiveVal("")
  
  # Current configuration (for preview and processing)
  # Note: No longer tracking "pending" vs "applied" - changes are immediate
  applied <- reactiveValues(
    use_iqr = FALSE,
    outlier_method = 'remove',
    outlier_detection = 'iqr',
    iqr_factors = character(0),
    dbscan_numeric = NULL,
    tech_agg = 'median'
  )
  
  # App launcher state
  launching_app <- reactiveVal(FALSE)
  
  # DBSCAN cache to avoid recalculating clusters
  dbscan_cache <- reactiveValues(
    data = NULL,
    eps = NULL,
    clusters = NULL
  )
  
  # Flag to prevent observers from firing during config loading
  loading_config <- reactiveVal(FALSE)
  
  # Helper: append to validation log
  append_validation <- function(...) {
    msg <- paste0(paste0(...), "\n")
    current <- validation_log()
    validation_log(paste0(current, msg))
  }
  
  return(list(
    # Core state
    validated = validated,
    generating_raw = generating_raw,
    validation_log = validation_log,
    applied = applied,
    
    # App launcher state
    launching_app = launching_app,
    
    # DBSCAN cache
    dbscan_cache = dbscan_cache,
    
    # Config loading flag
    loading_config = loading_config,
    
    # Helper functions
    append_validation = append_validation
  ))
}

