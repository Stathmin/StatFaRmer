# Spline utilities (time derivation, adaptive k, fit/backoff, emmeans at time)

library(tidyverse)
library(lme4)
library(emmeans)

#' Derive numeric time from available columns with validation and repair
#' @param data data.frame with time columns
#' @return data.frame with time_numeric column added (or unchanged if no time column found)
deriveTimeNumeric <- function(data) {
  source_used <- NA_character_
  if ('timestamp' %in% names(data)) {
    ts <- try(as.POSIXct(data$timestamp, tz = 'UTC'), silent = TRUE)
    if (!inherits(ts, 'try-error') && any(!is.na(ts))) {
      t0 <- suppressWarnings(min(ts, na.rm = TRUE))
      data$time_numeric <- as.numeric(difftime(ts, t0, units = 'days'))
      source_used <- 'timestamp'
    }
  }
  if (!('time_numeric' %in% names(data)) && 'timestamp_group' %in% names(data)) {
    tg <- try(as.POSIXct(data$timestamp_group, tz = 'UTC'), silent = TRUE)
    if (!inherits(tg, 'try-error') && any(!is.na(tg))) {
      t0 <- suppressWarnings(min(tg, na.rm = TRUE))
      data$time_numeric <- as.numeric(difftime(tg, t0, units = 'days'))
      source_used <- 'timestamp_group'
    }
  }
  if (!('time_numeric' %in% names(data)) && 'dbscan_cluster' %in% names(data)) {
    if (is.factor(data$dbscan_cluster)) {
      idx <- as.integer(data$dbscan_cluster)
    } else {
      idx <- match(data$dbscan_cluster, sort(unique(data$dbscan_cluster)))
    }
    data$time_numeric <- as.numeric(idx)
    source_used <- 'dbscan_cluster_order'
  }
  if ('unit' %in% names(data) && 'time_numeric' %in% names(data)) {
    data <- data %>% dplyr::group_by(.data$unit) %>% dplyr::mutate(
      time_numeric = rank(.data$time_numeric, ties.method = 'average')
    ) %>% dplyr::ungroup()
  }
  if (exists('logInfo', mode = 'function')) {
    logInfo(paste0('deriveTimeNumeric: source=', source_used))
  }
  data  # Return data directly, not list
}

#' Choose spline degrees of freedom k adaptively
#' @param n_time integer number of unique time points
#' @param k_min minimum k (default 3)
#' @param k_max maximum k (default 6)
#' @return integer k clamped to [k_min, k_max]
#' @details Rule: k = clamp(n_time, k_min, k_max)
chooseSplineK <- function(n_time, k_min = 3L, k_max = 6L) {
  if (is.null(n_time) || is.na(n_time) || n_time <= 0) return(k_min)
  # Clamp n_time to bounds: allows full flexibility up to k_max
  k <- max(k_min, min(as.integer(n_time), k_max))
  as.integer(k)
}

#' Fit spline model with controlled backoff ladder
#' @param fixed_formula_str fixed-effects formula string
#' @param data data.frame containing time_numeric
#' @param random_terms character vector of random terms like '(1|unit)'
#' @param allow_random_slopes whether to try random slopes on time
#' @param k spline degrees of freedom
#' @return list(model=model, model_type='lmer'|'lm', formula_used=character)
fitSplineModel <- function(fixed_formula_str, data, random_terms, allow_random_slopes, k) {
  rs_term <- ''
  if (allow_random_slopes && 'unit' %in% names(data)) {
    enough <- try(any(table(data$unit) > (k + 2), na.rm = TRUE), silent = TRUE)
    if (!inherits(enough, 'try-error') && isTRUE(enough)) {
      rs_term <- paste0(' + (splines::bs(time_numeric, k=', k, ')|unit)')
    }
  }
  base_spline <- paste0('splines::bs(time_numeric, k=', k, ')')
  spline_formula_str <- paste0(
    fixed_formula_str,
    ' + ', base_spline,
    rs_term,
    if (length(random_terms) > 0) paste0(' + ', paste(random_terms, collapse = ' + ')) else ''
  )
  ctrl <- lme4::lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)
  model <- try(suppressWarnings(lme4::lmer(as.formula(spline_formula_str), data = data, REML = TRUE, control = ctrl)), silent = TRUE)
  model_type <- 'lmer'
  formula_used <- spline_formula_str
  if (inherits(model, 'try-error') || lme4::isSingular(model, tol = 1e-5)) {
    spline_formula_str2 <- paste0(
      fixed_formula_str, ' + ', base_spline,
      if (length(random_terms) > 0) paste0(' + ', paste(random_terms, collapse = ' + ')) else ''
    )
    model2 <- try(suppressWarnings(lme4::lmer(as.formula(spline_formula_str2), data = data, REML = TRUE, control = ctrl)), silent = TRUE)
    if (!inherits(model2, 'try-error') && !lme4::isSingular(model2, tol = 1e-5)) {
      model <- model2
      formula_used <- spline_formula_str2
    } else {
      lm_formula <- as.formula(paste0(fixed_formula_str, ' + ', base_spline))
      model <- stats::lm(lm_formula, data = data)
      model_type <- 'lm'
      formula_used <- deparse(lm_formula)
    }
  }
  list(model = model, model_type = model_type, formula_used = formula_used)
}

#' emmeans at representative time points helper
#' @param model fitted model (lmerMod or lm)
#' @param times character vector among c('early','mid','late') or numeric values
#' @param spec emmeans spec formula string like '~ treatment | cultivar'
#' @return data.frame with emmeans results
emmeansAt <- function(model, times = c('early', 'mid', 'late'), spec = '~ treatment | cultivar') {
  if (inherits(model, 'lmerMod')) {
    tn <- model@frame$time_numeric
  } else if (!is.null(model$model) && 'time_numeric' %in% names(model$model)) {
    tn <- model$model$time_numeric
  } else {
    stop('emmeansAt requires time_numeric in model frame')
  }
  qmap <- function(label) {
    switch(label,
      early = as.numeric(stats::quantile(tn, 0.1, na.rm = TRUE)),
      mid = as.numeric(stats::quantile(tn, 0.5, na.rm = TRUE)),
      late = as.numeric(stats::quantile(tn, 0.9, na.rm = TRUE)),
      suppressWarnings(as.numeric(label))
    )
  }
  pts <- vapply(times, qmap, numeric(1))
  pts <- unique(pts[!is.na(pts)])
  re_form <- if (inherits(model, 'lmerMod')) NA else NULL
  res_list <- lapply(pts, function(tval) {
    emmeans::emmeans(model, as.formula(spec), at = list(time_numeric = tval), re.form = re_form)
  })
  dfs <- lapply(seq_along(res_list), function(i) {
    df <- as.data.frame(res_list[[i]])
    df$time_numeric <- pts[i]
    df
  })
  dplyr::bind_rows(dfs)
}


