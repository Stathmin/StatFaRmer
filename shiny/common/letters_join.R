constructInteractionKey <- function(df, factorName) {
  if (!grepl(':', factorName, fixed = TRUE)) {
    return(as.character(df[[factorName]]))
  }
  comps <- strsplit(factorName, ':', fixed = TRUE)[[1]]
  comps <- comps[comps %in% names(df)]
  if (length(comps) == 0) return(character(nrow(df)))
  do.call(paste, c(lapply(comps, function(cn) as.character(df[[cn]])), list(sep = ':')))
}

joinLettersToData <- function(df, lettersDf, factorName, byNames = NULL) {
  if (is.null(lettersDf) || nrow(lettersDf) == 0) {
    return(list(df = df, matched = 0L, unmatched = nrow(df)))
  }
  # Build a stable join key without mutating original factor column
  key_col <- constructInteractionKey(df, factorName)
  temp_col <- paste0('_', factorName, '_letter')
  df[['__join_key__']] <- as.character(key_col)
  # Ensure lettersDf group is character for robust join
  lettersDf$group <- as.character(lettersDf$group)
  
  cat("DEBUG joinLettersToData: factorName =", factorName, "\n")
  cat("DEBUG joinLettersToData: df nrow =", nrow(df), "lettersDf nrow =", nrow(lettersDf), "\n")
  cat("DEBUG joinLettersToData: df key sample =", paste(head(unique(df[['__join_key__']]), 3), collapse=", "), "\n")
  cat("DEBUG joinLettersToData: lettersDf group sample =", paste(head(unique(lettersDf$group), 3), collapse=", "), "\n")
  cat("DEBUG joinLettersToData: byNames =", paste(byNames, collapse=", "), "\n")
  
  # Build join mapping, including by/stratum columns when provided
  by_map <- c("__join_key__" = 'group')
  if (!is.null(byNames) && length(byNames) > 0) {
    for (bn in byNames) {
      if (bn %in% names(df) && bn %in% names(lettersDf)) {
        by_map[[bn]] <- bn
        cat("DEBUG joinLettersToData: added by column", bn, "\n")
      }
    }
  }
  cat("DEBUG joinLettersToData: by_map =", paste(names(by_map), "=", by_map, collapse=", "), "\n")
  
  out <- dplyr::left_join(df, lettersDf, by = by_map)
  out[[temp_col]] <- out$letter
  out$letter <- NULL
  # Clean helper key
  out[['__join_key__']] <- NULL
  matched <- sum(!is.na(out[[temp_col]]))
  cat("DEBUG joinLettersToData: matched =", matched, "unmatched =", nrow(out) - matched, "\n")
  list(df = out, matched = matched, unmatched = nrow(out) - matched)
}


