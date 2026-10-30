# Collapse plant × time panels to one row per plant over selected times.
# Metric formulas align with Growth Summaries (trapezoidal AUC, peak max).

.sf_plant_id_cols <- function(data) {
  candidates <- c(
    "unit", "V.T.R", "cultivar", "treatment", "repeat",
    "T:X:Y", "table_no", "pos_x", "pos_y"
  )
  intersect(candidates, names(data))
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

#' Collapse filtered times to one row per plant with selected metrics.
#'
#' @param data Plant × time data.frame.
#' @param trait Name of numeric trait column.
#' @param time_col Name of time/cluster column used for selection.
#' @param time_values Values of `time_col` to keep.
#' @param metrics Character vector among `"mean"`, `"auc"`, `"peak"`.
#' @return One row per `unit` (or `V.T.R` if no `unit`) with metric + design cols.
collapsePlantsBySelectedTimes <- function(data, trait, time_col, time_values, metrics = c("mean", "auc", "peak")) {
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
  allowed <- c("mean", "auc", "peak")
  bad <- setdiff(metrics, allowed)
  if (length(bad) > 0L) {
    stop("Unknown metrics: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (length(metrics) < 1L) {
    stop("`metrics` must contain at least one of: mean, auc, peak", call. = FALSE)
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
    rows[[i]] <- out
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
