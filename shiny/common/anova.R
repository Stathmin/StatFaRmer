# ANOVA and Mixed Models (delegation only)

library(tidyverse)
library(broom)
library(lme4)

#' Gate ANOVA before fitting: unfit designs must not run aov/lmer.
#' @param data Filtered data.frame
#' @param factors Character vector of requested ANOVA factors
#' @param response Response column name
#' @return list(ok, reason, factors, cardinality)
assessAnovaFitness <- function(data, factors, response) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 2L) {
    return(list(ok = FALSE, reason = "Need at least 2 rows in filtered data.", factors = character(0), cardinality = integer(0)))
  }
  if (is.null(response) || !nzchar(response) || !response %in% names(data)) {
    return(list(ok = FALSE, reason = "Response trait missing from filtered data.", factors = character(0), cardinality = integer(0)))
  }
  y <- suppressWarnings(as.numeric(data[[response]]))
  if (sum(is.finite(y)) < 2L) {
    return(list(ok = FALSE, reason = "Response has fewer than 2 finite values.", factors = character(0), cardinality = integer(0)))
  }
  factors <- as.character(factors)
  factors <- factors[!is.na(factors) & nzchar(factors)]
  factors <- factors[factors %in% names(data)]
  if (length(factors) < 1L) {
    return(list(
      ok = FALSE,
      reason = "No ANOVA factors selected (or none exist in this project data). Intercept-only ANOVA is not run.",
      factors = character(0),
      cardinality = integer(0)
    ))
  }
  cardinality <- vapply(factors, function(f) {
    length(unique(stats::na.omit(data[[f]])))
  }, integer(1))
  usable <- factors[cardinality >= 2L]
  if (length(usable) < 1L) {
    detail <- paste(paste0(names(cardinality), "=", cardinality), collapse = ", ")
    return(list(
      ok = FALSE,
      reason = paste0(
        "ANOVA not run: every selected factor has <2 levels in the filtered data (",
        detail, "). Change filters/factors, or use BLUP / Growth Summaries."
      ),
      factors = character(0),
      cardinality = cardinality
    ))
  }
  list(ok = TRUE, reason = NULL, factors = usable, cardinality = cardinality)
}

#' Analyze cardinality of factors in data
#' @param data Filtered data
#' @param factors Character vector of factor names to check
#' @return Named integer vector of unique level counts
analyzeFactorCardinality <- function(data, factors) {
  cat("DEBUG analyzeFactorCardinality: factors =", paste(factors, collapse=", "), "\n")
  cat("DEBUG analyzeFactorCardinality: data columns =", paste(names(data), collapse=", "), "\n")
  result <- sapply(factors, function(f) {
    if (f %in% names(data)) {
      n_unique <- length(unique(data[[f]]))
      n_levels <- if (is.factor(data[[f]])) nlevels(data[[f]]) else NA
      cat("DEBUG analyzeFactorCardinality: factor", f, "has", n_unique, "unique values\n")
      if (!is.na(n_levels) && n_unique != n_levels) {
        cat(sprintf("DEBUG analyzeFactorCardinality: %s has %d unique values but %d factor levels\n", 
                    f, n_unique, n_levels))
      }
      n_unique
    } else {
      cat("DEBUG analyzeFactorCardinality: factor", f, "not found in data\n")
      0L
    }
  })
  cat("DEBUG analyzeFactorCardinality: result =", paste(names(result), "=", result, collapse=", "), "\n")
  result
}

#' Run classical ANOVA (aov) with given formula and data
#' @param data data.frame
#' @param formula character or formula
#' @return list(model=aov, anova_table=tibble)
runClassicANOVA <- function(data, formula) {
  f <- if (inherits(formula, 'formula')) formula else as.formula(formula)
  model <- aov(f, data = data)
  anova_table <- broom::tidy(model)
  list(model = model, anova_table = anova_table)
}

#' Placeholder delegator for performANOVA to be wired via model_resolver
#' Keeps backward compatibility until stats.R is refactored to call resolver
#' @param data data.frame
#' @param formula character/formula
#' @param out_variable response variable name
#' @param force_method optional override
#' @return list with model, model_type, anova_table, desc_stats, cardinality, blocking_factors
performANOVA <- function(data, formula, out_variable, force_method = NULL) {
  cat("DEBUG performANOVA: starting with formula =", formula, "\n")
  cat("DEBUG performANOVA: out_variable =", out_variable, "\n")
  cat("DEBUG performANOVA: data nrows =", nrow(data), "\n")
  
  # Temporary: classical ANOVA only; resolver wiring will extend this
  res <- runClassicANOVA(data, formula)
  cat("DEBUG performANOVA: runClassicANOVA completed\n")
  # Build minimal compatible structure; other fields left minimal for now
  cardinality <- tryCatch({
    formula_obj <- as.formula(formula)
    all_factors <- all.vars(formula_obj)[-1]
    analyzeFactorCardinality(data, all_factors)
  }, error = function(e) integer(0))
  desc_stats <- tryCatch({
    all_factors <- if (length(cardinality) > 0) names(cardinality) else character(0)
    cat("DEBUG desc_stats: all_factors =", paste(all_factors, collapse=", "), "\n")
    cat("DEBUG desc_stats: cardinality =", paste(names(cardinality), "=", cardinality, collapse=", "), "\n")
    if (length(all_factors) > 0) {
      result <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(all_factors))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean = mean(.data[[out_variable]], na.rm = TRUE),
          sd = sd(.data[[out_variable]], na.rm = TRUE),
          se = sd / sqrt(n),
          .groups = 'drop'
        )
      cat("DEBUG desc_stats: result rows =", nrow(result), "\n")
      result
    } else {
      cat("DEBUG desc_stats: no factors, returning empty tibble\n")
      tibble::tibble()
    }
  }, error = function(e) {
    cat("DEBUG desc_stats: ERROR =", e$message, "\n")
    tibble::tibble()
  })
  result <- list(
    model = res$model,
    model_type = 'aov',
    anova_table = res$anova_table,
    desc_stats = desc_stats,
    cardinality = cardinality,
    blocking_factors = character(0),
    success = TRUE
  )
  cat("DEBUG performANOVA: returning cardinality =", paste(names(cardinality), "=", cardinality, collapse=", "), "\n")
  result
}


