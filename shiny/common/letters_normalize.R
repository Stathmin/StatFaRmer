normalizeGroupLabels <- function(groups, factorName) {
  if (length(groups) == 0) return(groups)
  comps <- strsplit(factorName, ':', fixed = TRUE)[[1]]
  vapply(groups, function(g) {
    parts <- strsplit(as.character(g), ':', fixed = TRUE)[[1]]
    if (length(parts) != length(comps)) return(as.character(g))
    paste(mapply(function(part, comp) sub(paste0('^', comp), '', part), parts, comps, USE.NAMES = FALSE), collapse = ':')
  }, character(1))
}


