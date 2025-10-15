shouldFilterComparisons <- function(nComparisons) {
  is.finite(nComparisons) && !is.na(nComparisons) && nComparisons >= 1000
}

filterSimpleComparisons <- function(tukeyDf) {
  if (is.null(tukeyDf) || nrow(tukeyDf) == 0 || !('comparison' %in% names(tukeyDf))) return(tukeyDf)
  parts <- strsplit(tukeyDf$comparison, '-')
  keep_mask <- vapply(parts, function(p) {
    if (length(p) != 2) return(FALSE)
    lhs <- trimws(p[[1]]); rhs <- trimws(p[[2]])
    lhs_levels <- strsplit(lhs, ':', fixed = TRUE)[[1]]
    rhs_levels <- strsplit(rhs, ':', fixed = TRUE)[[1]]
    if (length(lhs_levels) != length(rhs_levels)) return(FALSE)
    sum(lhs_levels != rhs_levels) == 1
  }, logical(1))
  tukeyDf[keep_mask, , drop = FALSE]
}


