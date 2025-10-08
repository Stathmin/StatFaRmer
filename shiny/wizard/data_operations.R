# =============================================================================
# WIZARD DATA OPERATIONS
# =============================================================================
# Centralized data loading, transformation, and processing functions
# 
# Note: loadProjectDataWithFallback() is now in common/utils.R for reuse

#' Get available variables for project with logit transform support
#' @param proj Project name
#' @param generating_raw Reactive value indicating if raw data is being generated
#' @param validated Reactive value indicating if project is validated
#' @param use_logit Whether logit transform is enabled
#' @param logit_treat_inf Whether to treat infinite values in logit
#' @return List with factors and numerics
getAvailableVariables <- function(proj, generating_raw, validated, use_logit, logit_treat_inf) {
  if (isTRUE(generating_raw())) return(list(factors = character(0), numerics = character(0)))
  if (!nzchar(proj) || !validated()) return(list(factors = character(0), numerics = character(0)))
  
  # Load raw data using unified function
  merged <- loadProjectDataWithFallback(proj, prefer = 'raw', context = 'availableVariables')
  if (is.null(merged)) return(list(factors = character(0), numerics = character(0)))
  
  if (nrow(merged) == 0) return(list(factors = character(0), numerics = character(0)))
  
  # Apply logit transform if enabled
  merged <- applyLogitTransform(merged, isTRUE(use_logit), isTRUE(logit_treat_inf))
  
  # Classify variables
  factors <- names(merged)[sapply(merged, function(x) is.character(x) || is.factor(x))]
  nums <- names(merged)[sapply(merged, is.numeric)]
  
  # Add dbscan_cluster to factors if it exists
  if ('dbscan_cluster' %in% names(merged) && !'dbscan_cluster' %in% factors) {
    factors <- c(factors, 'dbscan_cluster')
  }
  
  return(list(factors = sort(factors), numerics = sort(nums)))
}

