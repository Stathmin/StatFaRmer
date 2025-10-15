buildEmmeansCld <- function(emm, factorName, byNames = NULL) {
  # Build CLD using emmeans::cld; output data.frame(group, letter[, by...])
  cld_tab <- tryCatch(emmeans::cld(emm, adjust = 'tukey'), error = function(e) NULL)
  if (is.null(cld_tab) || nrow(cld_tab) == 0 || !('.group' %in% names(cld_tab))) {
    return(data.frame(group = character(0), letter = character(0), stringsAsFactors = FALSE))
  }
  # Comparison levels come from factorName column in emm grid
  grp <- as.character(cld_tab[[factorName]])
  letters_clean <- gsub('[[:space:]]+', '', as.character(cld_tab$.group))
  out <- data.frame(group = grp, letter = letters_clean, stringsAsFactors = FALSE)
  # Include by/stratum columns if provided and present in cld_tab
  if (!is.null(byNames) && length(byNames) > 0) {
    for (bn in byNames) {
      if (bn %in% names(cld_tab)) {
        out[[bn]] <- as.character(cld_tab[[bn]])
      }
    }
  }
  out
}


