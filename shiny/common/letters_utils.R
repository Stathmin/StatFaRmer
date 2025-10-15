lettersAreBlank <- function(lettersDf) {
  if (is.null(lettersDf) || nrow(lettersDf) == 0) return(TRUE)
  if (!('letter' %in% names(lettersDf))) return(TRUE)
  all(is.na(lettersDf$letter) | trimws(lettersDf$letter) == '')
}


relabelLettersByMeans <- function(emm, factorName, byNames = NULL, lettersDf) {
  # Ensure inputs are valid
  if (is.null(lettersDf) || nrow(lettersDf) == 0) return(lettersDf)
  grid <- tryCatch(as.data.frame(emm), error = function(e) NULL)
  if (is.null(grid) || !(factorName %in% names(grid))) return(lettersDf)
  emmean_col <- if ('emmean' %in% names(grid)) 'emmean' else if ('estimate' %in% names(grid)) 'estimate' else NULL
  if (is.null(emmean_col)) return(lettersDf)
  strata_cols <- if (!is.null(byNames) && length(byNames) > 0) intersect(byNames, names(grid)) else character(0)
  # Split by strata for deterministic relabeling
  if (length(strata_cols) > 0) {
    grid_splits <- split(grid, grid[strata_cols], drop = TRUE)
  } else {
    grid_splits <- list(`__all__` = grid)
  }
  out_list <- list()
  for (nm in names(grid_splits)) {
    gsub <- grid_splits[[nm]]
    # Determine stratum filter for lettersDf
    if (length(strata_cols) > 0 && nm != '__all__') {
      # Derive values from gsub (unique per stratum)
      vals <- lapply(strata_cols, function(sc) unique(as.character(gsub[[sc]]))[1])
      filt <- rep(TRUE, nrow(lettersDf))
      for (i in seq_along(strata_cols)) {
        sc <- strata_cols[i]
        v <- vals[[i]]
        if (sc %in% names(lettersDf)) {
          filt <- filt & (as.character(lettersDf[[sc]]) == v)
        }
      }
      lsub <- lettersDf[filt, , drop = FALSE]
    } else {
      lsub <- lettersDf
    }
    if (nrow(lsub) == 0) next
    # Join lsub with gsub to get means per group
    join_cols <- factorName
    names(gsub)[names(gsub) == factorName] <- 'group'
    means_df <- gsub[, c('group', emmean_col), drop = FALSE]
    names(means_df)[2] <- 'emmean'
    lsub$group <- as.character(lsub$group)
    merged <- merge(lsub, means_df, by = 'group', all.x = TRUE, sort = FALSE)
    if (nrow(merged) == 0) { out_list[[length(out_list) + 1]] <- lsub; next }
    # Rank letter clusters by max emmean (descending), then assign 'a','b',...
    if (!('letter' %in% names(merged))) { out_list[[length(out_list) + 1]] <- lsub; next }
    # If all emmeans are NA, skip relabel
    if (all(is.na(merged$emmean))) { out_list[[length(out_list) + 1]] <- lsub; next }
    agg <- stats::aggregate(emmean ~ letter, data = merged[!is.na(merged$emmean), , drop = FALSE], FUN = function(x) max(x, na.rm = TRUE))
    if (nrow(agg) == 0) { out_list[[length(out_list) + 1]] <- lsub; next }
    agg <- agg[order(-agg$emmean), , drop = FALSE]
    new_levels <- letters[seq_len(nrow(agg))]
    names(new_levels) <- as.character(agg$letter)
    merged$letter <- as.character(new_levels[as.character(merged$letter)])
    # Restore factorName column name in gsub copy
    names(gsub)[names(gsub) == 'group'] <- factorName
    out_list[[length(out_list) + 1]] <- merged
  }
  # Combine and keep original column order
  if (length(out_list) == 0) return(lettersDf)
  out <- do.call(rbind, out_list)
  out
}


