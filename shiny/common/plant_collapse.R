# Collapse plant × time panels to one row per plant over selected times.
# Metric formulas align with Growth Summaries (trapezoidal AUC, peak max,
# growth_rate = Δy/Δt, last = endpoint at last selected time).
#
# These named metrics are the ONLY allowed across-time summaries.
# Never use a silent raw-trait mean for violins/Tukey/posthoc OR as an
# unlabeled multi-trait DR phenotype matrix — see
# .sf_forbid_across_time_raw_mean_posthoc / .sf_dr_assert_named_temporal_metrics.
# For genotype distributions of a raw phenotype, slice ONE time instead.

.sf_plant_id_cols <- function(data) {
  candidates <- c(
    "unit", "V.T.R", "cultivar", "treatment", "repeat",
    "T:X:Y", "table_no", "pos_x", "pos_y"
  )
  intersect(candidates, names(data))
}

#' Apply left-sidebar cultivar / treatment filters to a plant-level table.
#' Shared by Growth Summaries consumers (Trait DR, BLUP) so tab plots respect
#' the same filters as `filteredData()` / main distPlot.
#' @param ps Plant-level data.frame (or NULL)
#' @param combined_inputs Shiny reactiveValues / list with treatments, cultivars
#' @return Filtered data.frame (or NULL / empty)
.sf_filter_plants_by_sidebar <- function(ps, combined_inputs) {
  if (is.null(ps) || !is.data.frame(ps) || nrow(ps) < 1L) {
    return(ps)
  }
  if (is.null(combined_inputs)) {
    return(ps)
  }
  treatments <- tryCatch(combined_inputs$treatments, error = function(e) NULL)
  cultivars <- tryCatch(combined_inputs$cultivars, error = function(e) NULL)
  out <- ps
  if (length(treatments) > 0L && "treatment" %in% names(out)) {
    out <- out[as.character(out$treatment) %in% as.character(treatments), , drop = FALSE]
  }
  if (length(cultivars) > 0L && "cultivar" %in% names(out)) {
    out <- out[as.character(out$cultivar) %in% as.character(cultivars), , drop = FALSE]
  }
  if (nrow(out) > 0L) {
    return(droplevels(out))
  }
  out
}

#' Snapshot of sidebar filters that affect plant-level sidecars.
.sf_sidebar_plant_snapshot <- function(combined_inputs) {
  if (is.null(combined_inputs)) {
    return(list(cultivars = character(0), treatments = character(0)))
  }
  cvs <- tryCatch(combined_inputs$cultivars, error = function(e) character(0))
  trt <- tryCatch(combined_inputs$treatments, error = function(e) character(0))
  if (is.null(cvs)) cvs <- character(0)
  if (is.null(trt)) trt <- character(0)
  list(
    cultivars = sort(unique(as.character(cvs))),
    treatments = sort(unique(as.character(trt)))
  )
}

#' Subset BLUP genotype ranks to current sidebar cultivars (display / plots).
.sf_blup_filter_geno_table <- function(geno_tbl, combined_inputs) {
  if (is.null(geno_tbl) || !is.data.frame(geno_tbl) || nrow(geno_tbl) < 1L) {
    return(geno_tbl)
  }
  snap <- .sf_sidebar_plant_snapshot(combined_inputs)
  if (length(snap$cultivars) < 1L || !"level" %in% names(geno_tbl)) {
    return(geno_tbl)
  }
  geno_tbl[as.character(geno_tbl$level) %in% snap$cultivars, , drop = FALSE]
}

#' Subset treatment EMMeans to current sidebar treatments.
.sf_blup_filter_trt_table <- function(trt_tbl, combined_inputs) {
  if (is.null(trt_tbl) || !is.data.frame(trt_tbl) || nrow(trt_tbl) < 1L) {
    return(trt_tbl)
  }
  snap <- .sf_sidebar_plant_snapshot(combined_inputs)
  if (length(snap$treatments) < 1L || !"level" %in% names(trt_tbl)) {
    return(trt_tbl)
  }
  trt_tbl[as.character(trt_tbl$level) %in% snap$treatments, , drop = FALSE]
}

#' Subset Dunnett contrast rows whose genotype levels remain in sidebar cultivars.
.sf_blup_filter_dunnett <- function(dunnett_df, combined_inputs) {
  if (is.null(dunnett_df) || !is.data.frame(dunnett_df) || nrow(dunnett_df) < 1L) {
    return(dunnett_df)
  }
  snap <- .sf_sidebar_plant_snapshot(combined_inputs)
  if (length(snap$cultivars) < 1L) {
    return(dunnett_df)
  }
  lab_col <- intersect(c("contrast", "Contrasts", "level", "cultivar"), names(dunnett_df))[1]
  if (is.na(lab_col)) {
    return(dunnett_df)
  }
  keep <- vapply(as.character(dunnett_df[[lab_col]]), function(lab) {
    # contrast labels like "G2 - G1" — keep if any non-check token is still selected
    toks <- trimws(unlist(strsplit(lab, "\\s*-\\s*|\\s+")))
    toks <- toks[nzchar(toks) & toks != "-"]
    any(toks %in% snap$cultivars)
  }, logical(1))
  dunnett_df[keep, , drop = FALSE]
}

.sf_metric_mean <- function(values) {
  mean(values, na.rm = TRUE)
}

.sf_metric_peak <- function(values) {
  if (all(is.na(values))) {
    return(NA_real_)
  }
  max(values, na.rm = TRUE)
}

#' Trapezoidal AUC using POSIXct/numeric time (hours), matching Growth Summaries.
.sf_metric_auc <- function(values, times) {
  ok <- !is.na(values) & !is.na(times)
  values <- values[ok]
  times <- times[ok]
  if (length(values) < 2L) {
    return(0)
  }
  ord <- order(times)
  values <- values[ord]
  times <- as.numeric(times[ord])
  sum(diff(times) / 3600 * (utils::head(values, -1) + utils::tail(values, -1)) / 2, na.rm = TRUE)
}

#' Value at the last selected time (endpoint), ordered by `times`.
.sf_metric_last <- function(values, times) {
  ok <- !is.na(values) & !is.na(times)
  if (!any(ok)) {
    return(NA_real_)
  }
  values <- values[ok]
  times <- times[ok]
  ord <- order(as.numeric(times))
  as.numeric(values[[ord[[length(ord)]]]])
}

#' Net growth rate (last − first) / Δt; Δt in days for POSIXt, else numeric axis units.
.sf_metric_growth_rate <- function(values, times) {
  ok <- !is.na(values) & !is.na(times)
  values <- values[ok]
  times <- times[ok]
  if (length(values) < 2L) {
    return(NA_real_)
  }
  ord <- order(as.numeric(times))
  values <- values[ord]
  times <- times[ord]
  if (inherits(times, "POSIXt")) {
    dt <- as.numeric(difftime(times[[length(times)]], times[[1L]], units = "days"))
  } else {
    dt <- as.numeric(times[[length(times)]]) - as.numeric(times[[1L]])
  }
  if (!is.finite(dt) || abs(dt) < .Machine$double.eps) {
    return(NA_real_)
  }
  (as.numeric(values[[length(values)]]) - as.numeric(values[[1L]])) / dt
}

.sf_collapse_allowed_metrics <- c("mean", "auc", "peak", "growth_rate", "last")

#' TRUE if column names look like named temporal metrics (GS or trait__metric).
.sf_is_named_temporal_metric_col <- function(nm) {
  nm <- as.character(nm)
  grepl("__(auc|peak|growth_rate|last|mean)$", nm) |
    nm %in% .sf_collapse_allowed_metrics
}

#' Abort if trait columns are unlabeled raw phenotypes (silent across-time means).
#' @param trait_cols character
#' @param context label for error
.sf_dr_assert_named_temporal_metrics <- function(trait_cols, context = "Trait DR") {
  trait_cols <- as.character(trait_cols)
  if (length(trait_cols) < 1L) {
    stop(context, ": no trait columns.", call. = FALSE)
  }
  bad <- trait_cols[!.sf_is_named_temporal_metric_col(trait_cols)]
  if (length(bad) > 0L) {
    stop(
      context, ": multi-trait DR requires named temporal metrics ",
      "(e.g. trait__auc / trait__peak / trait__growth_rate / trait__last), ",
      "not unlabeled raw phenotypes that imply silent across-time means. Bad: ",
      paste(utils::head(bad, 8L), collapse = ", "),
      if (length(bad) > 8L) "…" else "",
      call. = FALSE
    )
  }
  invisible(trait_cols)
}

#' Collapse filtered times to one row per plant with selected metrics.
#'
#' @param data Plant × time data.frame.
#' @param trait Name of numeric trait column.
#' @param time_col Name of time/cluster column used for selection.
#' @param time_values Values of `time_col` to keep.
#' @param metrics Character vector among `"mean"`, `"auc"`, `"peak"`, `"growth_rate"`, `"last"`.
#' @return One row per `unit` (or `V.T.R` if no `unit`) with metric + design cols.
collapsePlantsBySelectedTimes <- function(data, trait, time_col, time_values,
                                          metrics = c("mean", "auc", "peak")) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  if (!trait %in% names(data)) {
    stop("Can't collapse plants: missing trait column `", trait, "`", call. = FALSE)
  }
  if (!time_col %in% names(data)) {
    stop("Can't collapse plants: missing time column `", time_col, "`", call. = FALSE)
  }
  metrics <- unique(as.character(metrics))
  allowed <- .sf_collapse_allowed_metrics
  bad <- setdiff(metrics, allowed)
  if (length(bad) > 0L) {
    stop("Unknown metrics: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (length(metrics) < 1L) {
    stop(
      "`metrics` must contain at least one of: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  # Soft-tolerate missing geometry: collapse itself does not require X/Y;
  # ensureTxyColumns fills them when parseable (T:X:Y / unit), else no-op.
  # Spatial residualization (not collapse) owns hard geometry requirements.
  data <- ensureTxyColumns(data)
  filtered <- data[data[[time_col]] %in% time_values, , drop = FALSE]
  if (nrow(filtered) < 1L) {
    stop("Can't collapse plants: no rows match selected time_values", call. = FALSE)
  }

  id_col <- if ("unit" %in% names(filtered)) {
    "unit"
  } else if ("V.T.R" %in% names(filtered)) {
    "V.T.R"
  } else {
    stop("Can't collapse plants: need `unit` or `V.T.R`", call. = FALSE)
  }

  design_cols <- .sf_plant_id_cols(filtered)
  time_for_auc <- if ("timestamp_group" %in% names(filtered)) {
    filtered$timestamp_group
  } else {
    filtered[[time_col]]
  }

  ids <- unique(filtered[[id_col]])
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    idx <- which(filtered[[id_col]] == ids[[i]])
    plant <- filtered[idx, , drop = FALSE]
    values <- plant[[trait]]
    times <- time_for_auc[idx]

    out <- plant[1L, design_cols, drop = FALSE]
    if ("mean" %in% metrics) {
      out$mean <- .sf_metric_mean(values)
    }
    if ("auc" %in% metrics) {
      out$auc <- .sf_metric_auc(values, times)
    }
    if ("peak" %in% metrics) {
      out$peak <- .sf_metric_peak(values)
    }
    if ("growth_rate" %in% metrics) {
      out$growth_rate <- .sf_metric_growth_rate(values, times)
    }
    if ("last" %in% metrics) {
      out$last <- .sf_metric_last(values, times)
    }
    rows[[i]] <- out
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Collapse many traits to one plant-level wide table of named temporal metrics.
#'
#' Columns are `trait__metric` (e.g. `digital_biomass_mm3__auc`). Suitable as
#' Trait DR input — never unlabeled across-time raw means.
#'
#' @param data Plant × time data.frame
#' @param traits Character vector of numeric trait columns
#' @param time_col Time/cluster column
#' @param time_values Levels to keep
#' @param metrics Among mean/auc/peak/growth_rate/last (DR default: no silent mean-only)
#' @return Wide plant data.frame + attr `metric_cols`
collapseManyTraitsBySelectedTimes <- function(data,
                                              traits,
                                              time_col,
                                              time_values,
                                              metrics = c("auc", "peak", "growth_rate", "last")) {
  traits <- unique(as.character(traits))
  traits <- intersect(traits, names(data))
  if (length(traits) < 1L) {
    stop("collapseManyTraitsBySelectedTimes: no traits found in data", call. = FALSE)
  }
  metrics <- unique(as.character(metrics))
  id_col <- if ("unit" %in% names(data)) {
    "unit"
  } else if ("V.T.R" %in% names(data)) {
    "V.T.R"
  } else {
    stop("collapseManyTraitsBySelectedTimes: need unit or V.T.R", call. = FALSE)
  }

  out <- NULL
  for (tr in traits) {
    one <- collapsePlantsBySelectedTimes(
      data = data,
      trait = tr,
      time_col = time_col,
      time_values = time_values,
      metrics = metrics
    )
    ren <- intersect(metrics, names(one))
    for (m in ren) {
      names(one)[names(one) == m] <- paste0(tr, "__", m)
    }
    keep <- c(intersect(.sf_plant_id_cols(one), names(one)), paste0(tr, "__", ren))
    one <- one[, keep, drop = FALSE]
    out <- if (is.null(out)) {
      one
    } else {
      merge(out, one, by = intersect(names(out), names(one)), all = TRUE)
    }
  }
  metric_cols <- grep("__", names(out), value = TRUE)
  .sf_dr_assert_named_temporal_metrics(metric_cols, context = "collapseManyTraitsBySelectedTimes")
  attr(out, "metric_cols") <- metric_cols
  attr(out, "traits") <- traits
  attr(out, "metrics") <- metrics
  attr(out, "time_col") <- time_col
  attr(out, "time_values") <- time_values
  out
}
