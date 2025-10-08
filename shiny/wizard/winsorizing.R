# =============================================================================
# WINSORIZING FUNCTIONS (CONSOLIDATED)
# =============================================================================
# Uses computeOutlierBounds() from outlier_operations.R for bounds calculation

#' Winsorize data (consolidated IQR and Z-score methods)
#' @param data Data frame
#' @param method Winsorizing method ('iqr' or 'zscore')
#' @param variables Variables to winsorize
#' @param factors Factor columns for grouping
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @return Data frame with winsorized values
winsorizeData <- function(data, method, variables, factors = character(0), 
                         iqr_multiplier = 1.5, zscore_threshold = 2.5) {
  if (length(variables) == 0) return(data)
  if (!method %in% c('iqr', 'zscore')) return(data)
  
  result <- data
  
  for (var in variables) {
    if (!var %in% names(data)) next
    
    if (length(factors) == 0) {
      # No grouping - winsorize across all data
      values <- data[[var]]
      bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
      
      if (!is.null(bounds)) {
        result[[var]] <- pmax(pmin(values, bounds$upper_bound), bounds$lower_bound)
      }
    } else {
      # Group by factors
      result <- result %>%
        dplyr::group_by(!!!rlang::syms(factors[factors %in% names(data)])) %>%
        dplyr::mutate(!!var := {
          values <- .data[[var]]
          bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
          
          if (!is.null(bounds)) {
            pmax(pmin(values, bounds$upper_bound), bounds$lower_bound)
          } else {
            values
          }
        }) %>%
        dplyr::ungroup()
    }
  }
  
  return(result)
}
