# Model resolver: decides between aov, lmer, and spline LMM

library(tidyverse)
library(broom)
library(car)
library(lme4)

# Guardrail defaults (can be overridden upstream if sourced before)
if (!exists('MAX_FIT_SECONDS')) MAX_FIT_SECONDS <- 5
if (!exists('MAX_ROWS_FOR_SPLINE')) MAX_ROWS_FOR_SPLINE <- 50000

#' Resolve and fit appropriate model (aov/lmer/spline)
#' @param data data.frame
#' @param formula character/formula
#' @param out_variable response variable name
#' @param force_method optional: 'aov','lmer','spline'
#' @return list with model, model_type, anova_table, desc_stats, cardinality, blocking_factors
resolveModel <- function(data, formula, out_variable, force_method = NULL) {
  t_start <- proc.time()[['elapsed']]
  # Prepare time
  if (exists('deriveTimeNumeric', mode = 'function')) {
    data <- deriveTimeNumeric(data)
  }
  # Parse factors
  formula_obj <- as.formula(formula)
  all_factors <- all.vars(formula_obj)[-1]
  # Ensure factor levels reflect current data subset for accurate cardinality
  data <- droplevels(data)
  cardinality <- analyzeFactorCardinality(data, all_factors)
  cat("DEBUG model_resolver: cardinality =", paste(names(cardinality), "=", cardinality, collapse=", "), "\n")
  # Hard override: force aov
  if (!is.null(force_method) && identical(force_method, 'aov')) {
    model <- aov(as.formula(formula), data = data)
    anova_table <- broom::tidy(model)
    desc_stats <- tryCatch({
      data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(all_factors))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean = mean(.data[[out_variable]], na.rm = TRUE),
          sd = sd(.data[[out_variable]], na.rm = TRUE),
          se = sd / sqrt(n),
          .groups = 'drop'
        )
    }, error = function(e) tibble::tibble())
    if (exists('logInfo', mode = 'function')) logInfo('resolveModel: forced aov override')
    return(list(
      model = model,
      model_type = 'aov',
      anova_table = anova_table,
      desc_stats = desc_stats,
      cardinality = cardinality,
      blocking_factors = character(0),
      success = TRUE
    ))
  }
  # Only time-like factors trigger mixed models; do NOT use 'unit' to trigger
  high_card_factors <- names(cardinality[cardinality > 10])
  time_like <- intersect(high_card_factors, c('dbscan_cluster', 'timestamp_group'))
  blocking_candidates <- time_like
  # Explicit fast-path: if dbscan_cluster exists and levels <= 10, prefer aov unless overridden
  if (is.null(force_method) && 'dbscan_cluster' %in% names(cardinality) && !is.na(cardinality[['dbscan_cluster']]) && cardinality[['dbscan_cluster']] <= 10) {
    model <- aov(as.formula(formula), data = data)
    anova_table <- broom::tidy(model)
    desc_stats <- tryCatch({
      data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(all_factors))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean = mean(.data[[out_variable]], na.rm = TRUE),
          sd = sd(.data[[out_variable]], na.rm = TRUE),
          se = sd / sqrt(n),
          .groups = 'drop'
        )
    }, error = function(e) tibble::tibble())
    if (exists('logInfo', mode = 'function')) logInfo('resolveModel: fast-path aov for <=10 dbscan_cluster levels')
    return(list(
      model = model,
      model_type = 'aov',
      anova_table = anova_table,
      desc_stats = desc_stats,
      cardinality = cardinality,
      blocking_factors = character(0),
      success = TRUE
    ))
  }
  use_spline <- (!is.null(force_method) && identical(force_method, 'spline')) ||
                ('dbscan_cluster' %in% names(cardinality) && cardinality[['dbscan_cluster']] > 10)
  use_lmer <- if (!is.null(force_method)) {
    force_method %in% c('lmer', 'spline')
  } else {
    length(blocking_candidates) > 0
  }
  # Hard skip for very high cardinality in auto mode: avoid attempting any lmer/spline fit
  max_cardinality <- if (length(cardinality) > 0) max(cardinality, na.rm = TRUE) else 0
  if (is.null(force_method) && is.finite(max_cardinality) && max_cardinality > 20) {
    use_lmer <- FALSE
    use_spline <- FALSE
    blocking_candidates <- character(0)
  }
  model_type <- if (use_lmer) 'lmer' else 'aov'
  # Fixed effects
  fixed_factors <- if (use_lmer) setdiff(all_factors, blocking_candidates) else all_factors
  if (length(fixed_factors) == 0) {
    fixed_formula_str <- paste(out_variable, '~ 1')
  } else if (length(fixed_factors) == 1) {
    fixed_formula_str <- paste(out_variable, '~', fixed_factors[1])
  } else {
    fixed_formula_str <- paste(out_variable, '~', paste0('(', paste(fixed_factors, collapse = ' + '), ')^2'))
  }
  # Random effects: only time-like blocking as random intercepts (unit is never used)
  random_str <- c()
  # Validate blocking terms; if none valid and not forced, disable lmer
  if (exists('buildValidatedRandomEffects', mode = 'function')) {
    re_val <- buildValidatedRandomEffects(data, blocking_candidates, include_unit = FALSE)
    if (!is.null(re_val$formula) && nzchar(re_val$formula)) {
      random_str <- unlist(strsplit(re_val$formula, ' \\+ '))
    } else {
      if (is.null(force_method)) {
        use_lmer <- FALSE
      }
      blocking_candidates <- character(0)
    }
  } else {
    for (bf in blocking_candidates) {
      if (bf %in% names(data)) random_str <- c(random_str, paste0('(1|', bf, ')'))
    }
  }
  # If lmer is forced but there are no random terms, create a small synthetic grouping factor
  if (!is.null(force_method) && identical(force_method, 'lmer') && length(random_str) == 0) {
    k_groups <- max(2L, min(5L, nrow(data) - 1L))
    data$.batch <- factor(rep(seq_len(k_groups), length.out = nrow(data)))
    random_str <- c(random_str, '(1|.batch)')
  }
  # Log selection inputs
  if (exists('logInfo', mode = 'function')) {
    logInfo(paste0('resolveModel: factors=', paste(all_factors, collapse=','),
                   ', cardinality=', paste(names(cardinality), cardinality, sep=':', collapse=','),
                   ', force_method=', ifelse(is.null(force_method),'NULL', force_method),
                   ', use_lmer=', use_lmer, ', use_spline=', use_spline))
  }

  # Fit
  if (use_lmer) {
    if (use_spline && 'time_numeric' %in% names(data)) {
      # Size guardrail
      if (nrow(data) > MAX_ROWS_FOR_SPLINE) {
        if (exists('logWarn', mode = 'function')) logWarn(paste0('resolveModel: dataset too large for spline (n=', nrow(data), '), falling back to aov'))
        model <- aov(as.formula(formula), data = data)
      } else {
        # Time guardrail via setTimeLimit - wrap entire spline path
        k <- if (exists('chooseSplineK', mode='function')) chooseSplineK(length(unique(data$time_numeric))) else 4L
        allow_rs <- 'unit' %in% names(data)
        
        setTimeLimit(elapsed = MAX_FIT_SECONDS, transient = TRUE)
        on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
        
        fit <- try({
          fitSplineModel(fixed_formula_str, data, random_str, allow_rs, k)
        }, silent = TRUE)
        
        # Reset time limit immediately after attempt
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        
        if (inherits(fit, 'try-error')) {
          err_msg <- as.character(fit)
          if (grepl('reached elapsed time limit|reached CPU time limit', err_msg)) {
            if (exists('logWarn', mode = 'function')) {
              logWarn(paste0('resolveModel: spline fit timeout (>', MAX_FIT_SECONDS, 's), falling back to aov'))
            }
          } else {
            if (exists('logWarn', mode = 'function')) {
              logWarn(paste0('resolveModel: spline fit failed: ', err_msg))
            }
          }
          model <- aov(as.formula(formula), data = data)
        } else {
          model <- fit$model
          model_type <- fit$model_type
          if (exists('logInfo', mode = 'function')) {
            logInfo(paste0('resolveModel: spline fit ok, model_type=', model_type, ', k=', k))
          }
        }
      }
    } else {
      lmer_formula_str <- if (length(random_str) > 0) paste(fixed_formula_str, '+', paste(random_str, collapse = ' + ')) else fixed_formula_str
      ctrl <- lme4::lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)
      model <- lme4::lmer(as.formula(lmer_formula_str), data = data, REML = TRUE, control = ctrl)
      model_type <- 'lmer'
      if (lme4::isSingular(model, tol = 1e-5)) {
        # Keep singular lmer to preserve mixed model semantics (no unit term to drop)
        if (exists('logWarn', mode = 'function')) logWarn('resolveModel: lmer singular; returning singular lmer model')
      }
    }
  } else {
    model <- aov(as.formula(formula), data = data)
    model_type <- 'aov'
  }
  # ANOVA table and final model_type assignment
  if (inherits(model, 'lmerMod') || inherits(model, 'merMod')) {
    model_type <- 'lmer'
    # For high-cardinality factors (>20 levels), car::Anova is too slow - skip to aov
    # BUT respect forced spline/lmer override (user explicitly requested it)
    max_cardinality <- if (length(cardinality) > 0) max(cardinality, na.rm = TRUE) else 0
    if (max_cardinality > 20 && !(!is.null(force_method) && force_method %in% c('lmer', 'spline'))) {
      if (exists('logWarn', mode = 'function')) {
        logWarn(paste0('resolveModel: high cardinality (', max_cardinality, 
                      ' levels), skipping car::Anova on lmer, using aov instead'))
      }
      # Refit as aov for fast ANOVA table
      model <- aov(as.formula(formula), data = data)
      model_type <- 'aov'
      anova_table <- broom::tidy(model)
      blocking_candidates <- character(0)  # Clear blocking factors since we fell back to aov
    } else {
      # Normal lmer path for reasonable cardinality OR forced override
      if (max_cardinality > 20 && exists('logInfo', mode = 'function')) {
        logInfo(paste0('resolveModel: high cardinality (', max_cardinality, 
                      ') but force_method=', force_method, ', keeping lmer/spline model'))
      }
      # Wrap car::Anova in time limit to prevent hangs
      setTimeLimit(elapsed = MAX_FIT_SECONDS * 2, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      anova_table_result <- try(broom::tidy(car::Anova(model, type = 'II', test.statistic = 'F')), silent = TRUE)
      if (inherits(anova_table_result, 'try-error')) {
        if (exists('logWarn', mode = 'function')) logWarn('resolveModel: car::Anova timeout, falling back to aov')
        model <- aov(as.formula(formula), data = data)
        model_type <- 'aov'
        anova_table <- broom::tidy(model)
        blocking_candidates <- character(0)  # Clear blocking factors since we fell back to aov
      } else {
        anova_table <- anova_table_result
      }
    }
  } else if (inherits(model, 'lm')) {
    model_type <- 'lm'
    anova_table <- broom::tidy(anova(model))
  } else {
    model_type <- 'aov'
    anova_table <- broom::tidy(model)
  }
  # Descriptive stats
  desc_stats <- tryCatch({
    cat("DEBUG model_resolver desc_stats: all_factors =", paste(all_factors, collapse=", "), "\n")
    cat("DEBUG model_resolver desc_stats: out_variable =", out_variable, "\n")
    cat("DEBUG model_resolver desc_stats: data nrows =", nrow(data), "\n")
    result <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(all_factors))) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean = mean(.data[[out_variable]], na.rm = TRUE),
        sd = sd(.data[[out_variable]], na.rm = TRUE),
        se = sd / sqrt(n),
        .groups = 'drop'
      )
    cat("DEBUG model_resolver desc_stats: result rows =", nrow(result), "\n")
    result
  }, error = function(e) {
    cat("DEBUG model_resolver desc_stats: ERROR =", e$message, "\n")
    tibble::tibble()
  })
  elapsed <- proc.time()[['elapsed']] - t_start
  if (exists('logInfo', mode = 'function')) {
    logInfo(paste0('resolveModel: completed model_type=', model_type, ', elapsed_s=', round(elapsed,3)))
  }
  list(
    model = model,
    model_type = model_type,
    anova_table = anova_table,
    desc_stats = desc_stats,
    cardinality = cardinality,
    blocking_factors = blocking_candidates,
    success = TRUE
  )
}


