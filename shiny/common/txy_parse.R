# Parse T:X:Y (or unit holding table:x:y) into table_no / pos_x / pos_y.

#' True when values look like `table:x:y` (three integer fields).
.sf_looks_like_txy <- function(x) {
  x <- as.character(x)
  sample <- x[!is.na(x) & nzchar(x)]
  if (length(sample) < 1L) {
    return(FALSE)
  }
  parts <- strsplit(utils::head(sample, 20L), ":", fixed = TRUE)
  all(vapply(parts, function(p) {
    length(p) == 3L && !anyNA(suppressWarnings(as.integer(p)))
  }, logical(1)))
}

#' Resolve a column that holds `table:x:y` strings.
#'
#' Prefer literal `T:X:Y` / `table:x:y`; else `unit` when it matches the pattern
#' (soy-325 merged tables store geometry in `unit` only).
.sf_txy_source_col <- function(data) {
  for (nm in c("T:X:Y", "table:x:y")) {
    if (nm %in% names(data) && .sf_looks_like_txy(data[[nm]])) {
      return(nm)
    }
  }
  if ("unit" %in% names(data) && .sf_looks_like_txy(data$unit)) {
    return("unit")
  }
  NULL
}

#' Ensure `table_no`, `pos_x`, and `pos_y` columns exist when parseable.
#'
#' If all three geometry columns are already present, returns `data` unchanged
#' (does not overwrite existing X/Y). Otherwise parses a `table:x:y` string
#' column (`T:X:Y`, `table:x:y`, or `unit` when it holds that pattern) and fills
#' only missing geometry columns. When no parseable source exists, returns
#' `data` unchanged (soft-skip) so non-spatial collapse/BLUP can proceed;
#' spatial residualization remains responsible for requiring X/Y.
#'
#' @param data A data.frame (or tibble) with plant rows.
#' @return `data`, possibly with integer `table_no`, `pos_x`, `pos_y` added.
ensureTxyColumns <- function(data) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  needed <- c("table_no", "pos_x", "pos_y")
  if (all(needed %in% names(data))) {
    return(data)
  }

  txy_col <- .sf_txy_source_col(data)
  if (is.null(txy_col)) {
    # Soft-skip: incomplete geometry, nothing to parse (non-spatial path OK).
    return(data)
  }

  parts <- strsplit(as.character(data[[txy_col]]), ":", fixed = TRUE)
  bad <- vapply(parts, function(p) length(p) != 3L, logical(1))
  if (any(bad)) {
    stop(
      "Can't parse `", txy_col, "`: expected table:x:y for every row",
      call. = FALSE
    )
  }
  mat <- do.call(rbind, lapply(parts, as.integer))
  if (!"table_no" %in% names(data)) {
    data$table_no <- mat[, 1L]
  }
  if (!"pos_x" %in% names(data)) {
    data$pos_x <- mat[, 2L]
  }
  if (!"pos_y" %in% names(data)) {
    data$pos_y <- mat[, 3L]
  }
  data
}
