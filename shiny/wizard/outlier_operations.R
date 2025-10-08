# =============================================================================
# WIZARD OUTLIER OPERATIONS
# =============================================================================
# Centralized outlier detection, counting, and handling functions

#' Apply outlier strategy (unified wrapper for all strategies)
#' 
#' Single function that handles all outlier strategies: count, remove, winsorize, keep.
#' Returns both the processed data and outlier statistics.
#' 
#' @param data Data frame
#' @param strategy 'remove', 'winsorize', or 'keep'
#' @param method Detection method ('iqr' or 'zscore')
#' @param variables Variables to process
#' @param factors Factor columns for grouping
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @return List with data, outlier_count, and outlier_summary
applyOutlierStrategy <- function(data, strategy = 'remove', method = 'iqr', 
                                 variables = character(0), factors = character(0),
                                 iqr_multiplier = 1.5, zscore_threshold = 2.5) {
  
  if (length(variables) == 0) {
    return(list(
      data = data,
      outlier_count = 0,
      outlier_summary = list()
    ))
  }
  
  # Count outliers before processing
  outlier_summary <- countOutliers(data, method, variables, factors, iqr_multiplier, zscore_threshold)
  outlier_count <- sum(unlist(outlier_summary), na.rm = TRUE)
  
  # Apply strategy
  if (strategy == 'winsorize') {
    data <- winsorizeData(data, method, variables, factors, iqr_multiplier, zscore_threshold)
    data$outlier <- FALSE  # No outliers after winsorizing
  } else if (strategy == 'remove') {
    data <- replaceOutliersWithNA(data, method, variables, factors, iqr_multiplier, zscore_threshold)
    data$outlier <- FALSE  # Per-cell removal, no row-level outliers
  } else {
    # 'keep' - do nothing
    data$outlier <- FALSE
  }
  
  return(list(
    data = data,
    outlier_count = outlier_count,
    outlier_summary = outlier_summary
  ))
}

#' Compute outlier bounds (SINGLE SOURCE OF TRUTH)
#' @param values Numeric vector
#' @param method Detection method ('iqr' or 'zscore')
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @return List with lower_bound and upper_bound, or NULL if cannot compute
computeOutlierBounds <- function(values, method, iqr_multiplier = 1.5, zscore_threshold = 2.5) {
  # Check if we have any valid values
  valid_values <- values[!is.na(values)]
  if (length(valid_values) == 0) {
    return(NULL)
  }
  
  if (method == 'iqr') {
    q1 <- quantile(values, 0.25, na.rm = TRUE)
    q3 <- quantile(values, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    return(list(
      lower_bound = q1 - iqr_multiplier * iqr,
      upper_bound = q3 + iqr_multiplier * iqr
    ))
  } else if (method == 'zscore') {
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    if (!is.na(sd_val) && sd_val > 0) {
      return(list(
        lower_bound = mean_val - zscore_threshold * sd_val,
        upper_bound = mean_val + zscore_threshold * sd_val
      ))
    }
  }
  return(NULL)
}

#' Count outliers in data (for reporting)
#' @param data Data frame
#' @param method Detection method ('iqr' or 'zscore')
#' @param variables Variables to analyze
#' @param factors Factor columns for grouping
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @return List with outlier counts by variable
countOutliers <- function(data, method, variables, factors = character(0), 
                         iqr_multiplier = 1.5, zscore_threshold = 2.5) {
  if (length(variables) == 0) return(list())
  
  outlier_summary <- list()
  
  for (var in variables) {
    if (!var %in% names(data)) next
    
    if (length(factors) == 0) {
      # No grouping - count outliers across all data
      values <- data[[var]]
      bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
      
      if (!is.null(bounds)) {
        var_outliers <- sum(values < bounds$lower_bound | values > bounds$upper_bound, na.rm = TRUE)
      } else {
        var_outliers <- 0
      }
    } else {
      # Group by factors - count outliers within groups
      var_outliers <- data %>%
        dplyr::group_by(!!!rlang::syms(factors[factors %in% names(data)])) %>%
        dplyr::summarise(
          outlier_count = {
            values <- .data[[var]]
            bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
            if (!is.null(bounds)) {
              sum(values < bounds$lower_bound | values > bounds$upper_bound, na.rm = TRUE)
            } else {
              0
            }
          },
          .groups = 'drop'
        ) %>%
        dplyr::pull(outlier_count) %>%
        sum(na.rm = TRUE)
    }
    
    outlier_summary[[var]] <- var_outliers
  }
  
  return(outlier_summary)
}

#' Replace outliers with NA (per-cell removal)
#' @param data Data frame
#' @param method Detection method ('iqr' or 'zscore')
#' @param variables Variables to process
#' @param factors Factor columns for grouping
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @return Data frame with outliers replaced by NA
replaceOutliersWithNA <- function(data, method, variables, factors = character(0), 
                                 iqr_multiplier = 1.5, zscore_threshold = 2.5) {
  if (length(variables) == 0) return(data)
  
  result <- data
  
  for (var in variables) {
    if (!var %in% names(data)) next
    
    if (length(factors) == 0) {
      # No grouping - process across all data
      values <- data[[var]]
      bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
      
      if (!is.null(bounds)) {
        outlier_mask <- values < bounds$lower_bound | values > bounds$upper_bound
        result[[var]][outlier_mask] <- NA
      }
      # If method is 'keep' or bounds are NULL, result[[var]] remains unchanged
    } else {
      # Group by factors
      result <- result %>%
        dplyr::group_by(!!!rlang::syms(factors[factors %in% names(data)])) %>%
        dplyr::mutate(!!var := {
          values <- .data[[var]]
          bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
          
          if (!is.null(bounds)) {
            outlier_mask <- values < bounds$lower_bound | values > bounds$upper_bound
            values[outlier_mask] <- NA
          }
          values
        }) %>%
        dplyr::ungroup()
    }
  }
  
  return(result)
}

#' Compute outlier flags for data (for plotting)
#' @param data Data frame
#' @param variable Variable name to check for outliers
#' @param method Detection method ('iqr' or 'zscore')
#' @param factors Factor columns for grouping
#' @param iqr_multiplier IQR multiplier for IQR method (default 1.5)
#' @param zscore_threshold Z-score threshold for Z-score method (default 2.5)
#' @param add_outlier_col If TRUE, adds 'outlier' column (default FALSE for backward compatibility)
#' @return Data frame with outlier flags
computeOutlierFlags <- function(data, variable, method = 'iqr', factors = character(0), 
                               iqr_multiplier = 1.5, zscore_threshold = 2.5, add_outlier_col = FALSE) {
  if (!variable %in% names(data)) {
    data$.is_inlier_tmp <- TRUE
    if (add_outlier_col) {
      data$outlier <- FALSE
      data$.is_inlier_tmp <- NULL
    }
    return(data)
  }
  
  if (length(factors) > 0) {
    # Group by factors
    factors_valid <- factors[factors %in% names(data)]
    if (length(factors_valid) > 0) {
      data <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(factors_valid))) %>%
        dplyr::mutate(
          .is_inlier_tmp = {
            values <- .data[[variable]]
            bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
            if (!is.null(bounds)) {
              (values >= bounds$lower_bound) & (values <= bounds$upper_bound)
            } else {
              TRUE
            }
          }
        ) %>%
        dplyr::ungroup()
    } else {
      data$.is_inlier_tmp <- TRUE
    }
  } else {
    # No grouping - compute outliers across all data
    data <- data %>%
      dplyr::mutate(
        .is_inlier_tmp = {
          values <- .data[[variable]]
          bounds <- computeOutlierBounds(values, method, iqr_multiplier, zscore_threshold)
          if (!is.null(bounds)) {
            (values >= bounds$lower_bound) & (values <= bounds$upper_bound)
          } else {
            TRUE
          }
        }
      )
  }
  
  # Optionally add outlier column
  if (add_outlier_col) {
    data$outlier <- !data$.is_inlier_tmp
    data$.is_inlier_tmp <- NULL
  }
  
  return(data)
}

#' Generate outlier detection configuration as JSON data
#' @param outlier_results Outlier detection results
#' @return List with configuration data for JSON serialization
generateOutlierConfig <- function(outlier_results) {
  if (is.null(outlier_results)) {
    return(list())
  }
  
  return(list(
    outlier_detection_method = outlier_results$method,
    outlier_detection_factors = outlier_results$factors,
    outlier_detection_variables = outlier_results$variables,
    total_outliers_detected = outlier_results$total_outliers,
    outliers_by_variable = outlier_results$summary
  ))
}
