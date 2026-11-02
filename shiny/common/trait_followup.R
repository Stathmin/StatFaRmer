# Trait follow-up after DR/BLUP: growth curves + single-time Tukey/violins.
#
# HARD RULE — never average a raw phenotype across timepoints and treat that
# mean as a distributional / posthoc sample. StatFaRmer observes the same
# biological units repeatedly; dates are developmental states, not replicates.
# Allowed collapses are only explicit named metrics (GS mean/AUC/peak/rate)
# via collapsePlantsBySelectedTimes / Growth Summaries.

#' Abort if code tries the across-time raw-mean antipattern for posthoc/violins.
#' @param context short label for the caller
#' @return never returns (always stops)
.sf_forbid_across_time_raw_mean_posthoc <- function(context = "trait follow-up") {
  stop(
    "StatFaRmer antipattern blocked (", context, "): ",
    "do not average raw phenotype values across timepoints for violins/Tukey/ANOVA ",
    "distributions. Use a single time slice, or Growth Summaries named metrics ",
    "(mean/AUC/peak/rate) from collapsePlantsBySelectedTimes. ",
    "Temporal observations of the same plants are not interchangeable replicates.",
    call. = FALSE
  )
}

#' Pick time column used for slicing temporal panels.
.sf_followup_time_col <- function(data) {
  if ("dbscan_cluster" %in% names(data)) return("dbscan_cluster")
  if ("timestamp_group" %in% names(data)) return("timestamp_group")
  if ("time_numeric" %in% names(data)) return("time_numeric")
  NULL
}

#' Distinct time levels available for follow-up slicing (character).
.sf_followup_time_levels <- function(data) {
  tc <- .sf_followup_time_col(data)
  if (is.null(tc)) return(character(0))
  as.character(sort(unique(stats::na.omit(data[[tc]]))))
}

#' Keep exactly one time level. Refuses multi-time averaging.
#' @param data observation-level data.frame
#' @param time_value level to keep (character/numeric); NULL → last level (with note)
#' @return list(data, time_col, time_value, n_times_available, note)
.sf_slice_single_time <- function(data, time_value = NULL) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
    return(list(
      data = NULL, time_col = NULL, time_value = NULL,
      n_times_available = 0L, note = "No data"
    ))
  }
  tc <- .sf_followup_time_col(data)
  if (is.null(tc)) {
    return(list(
      data = data, time_col = NULL, time_value = NULL,
      n_times_available = 1L,
      note = "No time column — treating rows as a single cross-section"
    ))
  }
  levels <- .sf_followup_time_levels(data)
  n_av <- length(levels)
  if (n_av < 1L) {
    return(list(
      data = NULL, time_col = tc, time_value = NULL,
      n_times_available = 0L, note = "No time levels"
    ))
  }
  if (is.null(time_value) || !nzchar(as.character(time_value))) {
    time_value <- levels[[length(levels)]]
    note <- sprintf(
      "Single time slice = last level (%s=%s). Not an across-time mean.",
      tc, time_value
    )
  } else {
    time_value <- as.character(time_value)
    if (!time_value %in% levels) {
      return(list(
        data = NULL, time_col = tc, time_value = time_value,
        n_times_available = n_av,
        note = paste("Requested time not in data:", time_value)
      ))
    }
    note <- sprintf("Single time slice (%s=%s).", tc, time_value)
  }
  # Compare as character to avoid factor/numeric mismatch
  keep <- as.character(data[[tc]]) == time_value
  list(
    data = data[keep, , drop = FALSE],
    time_col = tc,
    time_value = time_value,
    n_times_available = n_av,
    note = note
  )
}

#' One row per plant at a single time (no across-time aggregation).
.sf_plants_at_time <- function(data, trait, group = "cultivar", time_value = NULL) {
  sl <- .sf_slice_single_time(data, time_value = time_value)
  if (is.null(sl$data) || nrow(sl$data) < 1L) {
    return(list(plant_df = NULL, slice = sl))
  }
  id <- if ("unit" %in% names(sl$data)) {
    "unit"
  } else if ("V.T.R" %in% names(sl$data)) {
    "V.T.R"
  } else {
    NULL
  }
  if (is.null(id) || !trait %in% names(sl$data) || !group %in% names(sl$data)) {
    return(list(plant_df = NULL, slice = sl))
  }
  df <- sl$data
  df[[trait]] <- suppressWarnings(as.numeric(df[[trait]]))
  df <- df[is.finite(df[[trait]]) & !is.na(df[[group]]), , drop = FALSE]
  if (nrow(df) < 1L) return(list(plant_df = NULL, slice = sl))

  # If duplicate rows per plant at same time, mean within that time only (sensor dups)
  keep <- unique(c(id, group, "treatment", "repeat", "table_no", trait))
  keep <- intersect(keep, names(df))
  plant <- df[, keep, drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(setdiff(keep, trait)))) |>
    dplyr::summarise(!!trait := mean(.data[[trait]], na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
  # Guard: caller must not have passed multi-time data without slicing
  if (!is.null(sl$time_col) && sl$n_times_available > 1L) {
    # sliced data should have 1 time level
    n_left <- length(unique(stats::na.omit(sl$data[[sl$time_col]])))
    if (n_left > 1L) {
      .sf_forbid_across_time_raw_mean_posthoc("plants_at_time: slice failed")
    }
  }
  list(plant_df = plant, slice = sl)
}

#' Intentionally removed: naive mean across times for posthoc.
#' @export
.sf_plant_trait_means <- function(...) {
  .sf_forbid_across_time_raw_mean_posthoc(".sf_plant_trait_means")
}

#' One-way ANOVA + Tukey CLD at a single timepoint (reuses performTukey).
#' @param time_value required when data has >1 time level (NULL = last level)
#' @return list(success, letters, tukey, reason, plant_df, slice)
computeTraitTukeyLetters <- function(data, trait, group = "cultivar",
                                     time_value = NULL) {
  got <- .sf_plants_at_time(data, trait, group = group, time_value = time_value)
  plant <- got$plant_df
  sl <- got$slice
  if (is.null(plant) || nrow(plant) < 3L) {
    return(list(
      success = FALSE, letters = NULL,
      reason = paste("Need ≥3 plants at one time.", sl$note %||% ""),
      plant_df = plant, slice = sl
    ))
  }
  nlev <- length(unique(stats::na.omit(plant[[group]])))
  if (nlev < 2L) {
    return(list(
      success = FALSE, letters = NULL,
      reason = "Need ≥2 group levels at this time.",
      plant_df = plant, slice = sl
    ))
  }
  # All-pairs Tukey / CLD is unusable for ~325 genotypes (~52k pairs).
  max_tukey_levels <- 20L
  if (nlev > max_tukey_levels) {
    return(list(
      success = FALSE, letters = NULL,
      reason = sprintf(
        "Tukey CLD skipped: %d genotypes > %d. Use BLUP ranks + Dunnett vs check; distribution plot shows top/bottom ends.",
        nlev, max_tukey_levels
      ),
      plant_df = plant, slice = sl
    ))
  }
  plant[[group]] <- factor(plant[[group]])
  frm <- stats::as.formula(paste0("`", trait, "` ~ `", group, "`"))
  fit <- tryCatch(stats::aov(frm, data = plant), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(
      success = FALSE, letters = NULL, reason = conditionMessage(fit),
      plant_df = plant, slice = sl
    ))
  }
  if (!exists("performTukey", mode = "function")) {
    return(list(
      success = FALSE, letters = NULL, reason = "performTukey not loaded.",
      plant_df = plant, slice = sl
    ))
  }
  tk <- tryCatch(performTukey(fit, group), error = function(e) e)
  if (inherits(tk, "error")) {
    return(list(
      success = FALSE, letters = NULL, reason = conditionMessage(tk),
      plant_df = plant, slice = sl
    ))
  }
  letters <- NULL
  if (isTRUE(tk$success) && !is.null(tk$results[[group]]$letters)) {
    letters <- tk$results[[group]]$letters
  }
  list(
    success = isTRUE(tk$success) && !is.null(letters) && nrow(letters) > 0L,
    letters = letters,
    tukey = tk,
    reason = if (isTRUE(tk$success)) NA_character_ else (tk$error %||% "Tukey failed"),
    plant_df = plant,
    slice = sl,
    model = fit
  )
}

#' Suggest follow-up traits from DR loadings (absolute PC1/PC2).
suggestFollowupTraits <- function(loadings_df, max_n = 6L) {
  if (is.null(loadings_df) || !is.data.frame(loadings_df) || !"trait" %in% names(loadings_df)) {
    return(character(0))
  }
  pcs <- intersect(c("PC1", "PC2"), names(loadings_df))
  if (length(pcs) < 1L) return(character(0))
  score <- abs(as.numeric(loadings_df[[pcs[[1]]]]))
  if (length(pcs) >= 2L) {
    score <- score + 0.5 * abs(as.numeric(loadings_df[[pcs[[2]]]]))
  }
  ord <- order(score, decreasing = TRUE)
  unique(as.character(loadings_df$trait[ord]))[seq_len(min(max_n, length(ord)))]
}

#' Build curve (full time) + violin(+Tukey) at a single time slice.
#' @param time_value time level for violin/Tukey; NULL = last available
buildTraitFollowupPlots <- function(obs_data,
                                    trait,
                                    group = "cultivar",
                                    time_value = NULL,
                                    run_tukey = TRUE) {
  if (is.null(obs_data) || !is.data.frame(obs_data) || !trait %in% names(obs_data)) {
    empty <- statfarmerEmptyPlot(paste("Cannot plot", trait))
    return(list(curve = empty, violin = empty, tukey = NULL, trait = trait))
  }
  tuk <- if (isTRUE(run_tukey)) {
    computeTraitTukeyLetters(obs_data, trait, group = group, time_value = time_value)
  } else {
    got <- .sf_plants_at_time(obs_data, trait, group = group, time_value = time_value)
    list(success = FALSE, letters = NULL, plant_df = got$plant_df, slice = got$slice)
  }
  letters <- if (!is.null(tuk) && isTRUE(tuk$success)) tuk$letters else NULL
  plant <- tuk$plant_df
  sl_note <- tuk$slice$note %||% ""
  tv <- tuk$slice$time_value %||% time_value

  list(
    trait = trait,
    time_value = tv,
    slice_note = sl_note,
    curve = createTraitGrowthCurvePlot(
      obs_data, trait, group = group,
      # Prefer cluster index over raw timestamps for multi-geno screens
      time_col = .sf_followup_time_col(obs_data)
    ),
    violin = createTraitViolinPlot(
      if (!is.null(plant)) plant else obs_data,
      trait,
      group = group,
      letters_df = letters,
      title = paste0(
        "Single-time distribution + Tukey: ", trait,
        if (!is.null(tv) && nzchar(as.character(tv))) paste0(" @ ", tv) else ""
      ),
      subtitle = paste0(
        sl_note,
        if (!is.null(letters)) " Tukey CLD letters." else ""
      )
    ),
    tukey = tuk
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
