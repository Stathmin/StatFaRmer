# Random Effects Validator
# Determines if random effects should be included based on data quality checks

## NOTE: unit is never modeled as a random effect in StatFaRmer.
## This function is retained for compatibility but always excludes unit.
shouldIncludeUnitRE <- function(data, unit_var = 'unit', min_obs_per_unit = 5, min_units = 3) {
  return(list(include = FALSE, reason = 'unit random effect disabled by design'))
}

#' Check if a blocking factor should be included as random effect
#' @param data data.frame
#' @param blocking_var character, name of blocking variable
#' @param min_levels integer, minimum number of levels (default 3)
#' @return list(include = TRUE/FALSE, reason = character)
shouldIncludeBlockingRE <- function(data, blocking_var, min_levels = 3) {
  
  # Check 1: Variable exists
  if (!blocking_var %in% names(data)) {
    return(list(include = FALSE, reason = paste(blocking_var, "not in data")))
  }
  
  # Check 2: Enough levels
  n_levels <- dplyr::n_distinct(data[[blocking_var]], na.rm = TRUE)
  if (n_levels < min_levels) {
    return(list(include = FALSE, reason = paste0("too few levels (", n_levels, " < ", min_levels, ")")))
  }
  
  # Check 3: Not all observations in one level (variance check)
  level_counts <- table(data[[blocking_var]])
  if (max(level_counts) / sum(level_counts) > 0.95) {
    return(list(include = FALSE, reason = "95%+ observations in single level"))
  }
  
  # Check 4: Reasonable level sizes (avoid levels with 1-2 observations)
  n_tiny_levels <- sum(level_counts < 3)
  if (n_tiny_levels / length(level_counts) > 0.5) {
    return(list(include = FALSE, reason = "50%+ levels have <3 observations"))
  }
  
  return(list(include = TRUE, reason = "passed all checks"))
}

#' Build random effects formula with validation
#' @param data data.frame
#' @param blocking_candidates character vector of potential blocking factors
#' @param include_unit logical, whether to consider unit RE
#' @return list(formula = character, included = character vector, excluded = list)
buildValidatedRandomEffects <- function(data, blocking_candidates, include_unit = TRUE) {
  
  random_terms <- c()
  included <- c()
  excluded <- list()
  
  # Check unit random effect
  # Unit RE disabled by design
  if (include_unit && 'unit' %in% names(data)) {
    excluded[['unit']] <- 'unit random effect disabled by design'
  }
  
  # Check blocking factors
  for (bf in blocking_candidates) {
    if (bf == 'unit') next  # Already handled
    
    block_check <- shouldIncludeBlockingRE(data, bf)
    if (block_check$include) {
      random_terms <- c(random_terms, paste0('(1|', bf, ')'))
      included <- c(included, bf)
    } else {
      excluded[[bf]] <- block_check$reason
    }
  }
  
  return(list(
    formula = if (length(random_terms) > 0) paste(random_terms, collapse = ' + ') else NULL,
    included = included,
    excluded = excluded
  ))
}

