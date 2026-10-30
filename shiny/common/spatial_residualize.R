# Per-table spatial residualization via mgcv (approach B stage 1).

#' Residualize a trait against per-table X/Y smooths.
#'
#' Fits `mgcv` smooths of the form T|X + T|Y
#' (`s(x, by = table) + s(y, by = table)`, coarse `k`). When `x_col` /
#' `y_pos_col` are absent, returns identity residuals (`y`) and
#' `meta$spatial_applied = FALSE`.
#'
#' @param data A data.frame with trait and (optionally) geometry columns.
#' @param y_col Name of the numeric trait column to residualize.
#' @param table_col Table / run block column (default `"table_no"`).
#' @param x_col Spatial X column (default `"pos_x"`).
#' @param y_pos_col Spatial Y column (default `"pos_y"`).
#' @param k Coarse basis dimension per by-smooth (default `4L`).
#' @param ... Passed to [mgcv::gam()] (e.g. `method`).
#' @return `list(data, meta)`. `data` gains `paste0(y_col, "_spat_resid")`.
#'   `meta` has `spatial_applied`, `formula`, `converged`, `per_table_n`,
#'   `resid_col`, `y_col`, and optional `message`.
residualizeSpatial <- function(data,
                               y_col,
                               table_col = "table_no",
                               x_col = "pos_x",
                               y_pos_col = "pos_y",
                               k = 4L,
                               ...) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  if (missing(y_col) || length(y_col) != 1L || !is.character(y_col)) {
    stop("y_col must be a single column name", call. = FALSE)
  }
  if (!y_col %in% names(data)) {
    stop("Can't residualize: missing trait column `", y_col, "`", call. = FALSE)
  }

  resid_col <- paste0(y_col, "_spat_resid")
  empty_n <- integer(0)
  names(empty_n) <- character(0)

  identity_result <- function(message = NULL) {
    out <- data
    out[[resid_col]] <- out[[y_col]]
    meta <- list(
      spatial_applied = FALSE,
      formula = NA_character_,
      converged = NA,
      per_table_n = empty_n,
      resid_col = resid_col,
      y_col = y_col,
      k = as.integer(k)
    )
    if (!is.null(message)) {
      meta$message <- message
    }
    list(data = out, meta = meta)
  }

  if (!x_col %in% names(data) || !y_pos_col %in% names(data)) {
    return(identity_result("X/Y columns absent; spatial stage skipped"))
  }
  if (!table_col %in% names(data)) {
    return(identity_result("table column absent; spatial stage skipped"))
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(identity_result("mgcv not available; spatial stage skipped"))
  }

  y <- data[[y_col]]
  if (!is.numeric(y)) {
    stop("Can't residualize: `", y_col, "` must be numeric", call. = FALSE)
  }

  model_df <- data.frame(
    .y = as.numeric(y),
    .x = as.numeric(data[[x_col]]),
    .yp = as.numeric(data[[y_pos_col]]),
    .table = data[[table_col]],
    stringsAsFactors = FALSE
  )
  ok <- stats::complete.cases(model_df)
  if (!any(ok)) {
    return(identity_result("no complete X/Y/table/y rows; spatial stage skipped"))
  }

  model_df$.table_f <- factor(model_df$.table)
  per_table_n <- as.integer(table(model_df$.table_f[ok]))
  names(per_table_n) <- levels(model_df$.table_f)

  # Coarse k: prefer requested k, but do not exceed unique positions when
  # a table has ≥3 levels; otherwise keep k≥3 and let mgcv shrink (e.g.
  # constant Y within a table on the mini fixture).
  ux <- vapply(
    split(model_df$.x[ok], model_df$.table_f[ok], drop = TRUE),
    function(z) length(unique(z)),
    integer(1)
  )
  uy <- vapply(
    split(model_df$.yp[ok], model_df$.table_f[ok], drop = TRUE),
    function(z) length(unique(z)),
    integer(1)
  )
  k_req <- as.integer(k)
  if (is.na(k_req) || k_req < 3L) {
    k_req <- 3L
  }
  clamp_k <- function(n_unique, k_req) {
    if (min(n_unique) >= 3L) {
      as.integer(min(k_req, min(n_unique)))
    } else {
      as.integer(min(k_req, max(min(n_unique), 3L)))
    }
  }
  k_x <- clamp_k(ux, k_req)
  k_y <- clamp_k(uy, k_req)

  fit_gam <- function(include_y) {
    if (include_y) {
      frm <- stats::as.formula(
        paste0(
          ".y ~ .table_f + s(.x, by = .table_f, k = ", k_x,
          ") + s(.yp, by = .table_f, k = ", k_y, ")"
        )
      )
    } else {
      frm <- stats::as.formula(
        paste0(".y ~ .table_f + s(.x, by = .table_f, k = ", k_x, ")")
      )
    }
    fit <- mgcv::gam(frm, data = model_df[ok, , drop = FALSE], ...)
    list(fit = fit, formula = frm)
  }

  fitted <- tryCatch(fit_gam(include_y = TRUE), error = function(e) e)
  if (inherits(fitted, "error")) {
    fitted <- tryCatch(fit_gam(include_y = FALSE), error = function(e) e)
  }
  if (inherits(fitted, "error")) {
    return(identity_result(paste0("mgcv fit failed: ", conditionMessage(fitted))))
  }

  fit <- fitted$fit
  resid_ok <- as.numeric(stats::residuals(fit))
  resid_full <- rep(NA_real_, nrow(data))
  resid_full[ok] <- resid_ok

  out <- data
  out[[resid_col]] <- resid_full

  converged <- NA
  if (!is.null(fit$converged)) {
    converged <- isTRUE(fit$converged)
  } else if (!is.null(fit$mgcv.conv)) {
    # mgcv stores convergence info in mgcv.conv (list or logical)
    conv <- fit$mgcv.conv
    if (is.logical(conv)) {
      converged <- isTRUE(conv)
    } else if (is.list(conv) && !is.null(conv$fully.converged)) {
      converged <- isTRUE(conv$fully.converged)
    } else {
      converged <- TRUE
    }
  } else {
    converged <- TRUE
  }

  formula_chr <- paste(deparse(fitted$formula, width.cutoff = 500L), collapse = " ")
  meta <- list(
    spatial_applied = TRUE,
    formula = formula_chr,
    converged = converged,
    per_table_n = per_table_n,
    resid_col = resid_col,
    y_col = y_col,
    k = k_req,
    k_x = as.integer(k_x),
    k_y = as.integer(k_y),
    n_used = as.integer(sum(ok)),
    n_total = as.integer(nrow(data))
  )
  list(data = out, meta = meta)
}
