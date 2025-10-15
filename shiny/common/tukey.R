# Tukey / emmeans utilities

library(emmeans)
library(multcompView)

#' Perform Tukey HSD / emmeans-based contrasts
#' @param anova_model ANOVA model object
#' @param tukey_factors Factors for Tukey test
#' @return List with Tukey results
performTukey <- function(anova_model, tukey_factors) {
  tryCatch({
    use_emmeans <- requireNamespace('emmeans', quietly = TRUE)
    # Preserve user order in the returned list
    tukey_results <- structure(list(), .Names = character(0))
    for (factor in tukey_factors) {
      if (use_emmeans) {
        is_lmer <- inherits(anova_model, 'lmerMod')
        if (is_lmer) {
          model_formula <- formula(anova_model)
          all_terms <- attr(terms(model_formula, fixed.only = TRUE), 'term.labels')
          if (is.null(all_terms) || length(all_terms) == 0) all_terms <- character(0)
          fixed_terms <- if (length(all_terms) > 0) all_terms[!grepl('\\|', all_terms)] else character(0)
          if (length(fixed_terms) > 0) {
            fixed_terms_char <- as.character(fixed_terms)
            all_factors_in_model <- unique(unlist(lapply(fixed_terms_char, function(x) strsplit(x, ':'))))
          } else {
            all_factors_in_model <- character(0)
          }
          if (!factor %in% all_factors_in_model) {
            next
          }
          other <- setdiff(intersect(tukey_factors, all_factors_in_model), factor)
        } else {
          trms <- attr(terms(anova_model), 'term.labels')
          other <- setdiff(trms, factor)
        }
        # For each Tukey factor, compute emmeans with other Tukey factors as strata (by)
        by_vars <- setdiff(tukey_factors, factor)
        # Build emmeans with better df method for lmer to avoid NA p-values
        if (length(by_vars) > 0) {
          if (is_lmer) {
            emm <- emmeans::emmeans(anova_model, specs = factor, by = by_vars, lmer.df = 'satterthwaite')
          } else {
            emm <- emmeans::emmeans(anova_model, specs = factor, by = by_vars)
          }
        } else {
          if (is_lmer) {
            emm <- emmeans::emmeans(anova_model, specs = factor, lmer.df = 'satterthwaite')
          } else {
            emm <- emmeans::emmeans(anova_model, specs = factor)
          }
        }
        pw <- emmeans::contrast(emm, method = 'pairwise', adjust = 'tukey')
        pw_df <- as.data.frame(pw)
        if (nrow(pw_df) > 0 && 'contrast' %in% names(pw_df)) {
          contrast_str <- as.character(pw_df$contrast)
          lhs_rhs <- strsplit(contrast_str, ' - ')
          lhs <- vapply(lhs_rhs, function(p) trimws(p[[1]]), character(1))
          rhs <- vapply(lhs_rhs, function(p) trimws(p[[2]]), character(1))
          if (identical(factor, 'dbscan_cluster')) {
            lhs <- paste('Cluster', lhs)
            rhs <- paste('Cluster', rhs)
          }
          pw_df$comparison <- paste0(factor, ': ', lhs, ' vs ', rhs)
        } else {
          pw_df$comparison <- character(0)
        }
        if (!'p.adj' %in% names(pw_df)) {
          if ('p.value' %in% names(pw_df)) pw_df$p.adj <- pw_df$p.value
        }
        # Build CLD via helpers; if blank, fallback to p-value based letters
        cld_df <- buildEmmeansCld(emm, factor, byNames = by_vars)
        if (is.null(cld_df) || nrow(cld_df) == 0 || lettersAreBlank(cld_df)) {
          cld_df <- buildLettersFromPvals(pw_df, factor, byNames = by_vars)
        }
        # Normalize letter ordering so that 'a' corresponds to highest emmean per stratum
        if (!is.null(cld_df) && nrow(cld_df) > 0 && !lettersAreBlank(cld_df)) {
          cld_df <- relabelLettersByMeans(emm, factor, by_vars, cld_df)
        }
        # Graceful degradation: only for strata with a single level of the factor
        if (is.null(cld_df) || nrow(cld_df) == 0 || lettersAreBlank(cld_df)) {
          grid_df <- tryCatch(as.data.frame(emm), error = function(e) NULL)
          if (!is.null(grid_df)) {
            assign_list <- list()
            if (length(by_vars) > 0) {
              splits <- split(grid_df, grid_df[by_vars], drop = TRUE)
            } else {
              splits <- list(`__all__` = grid_df)
            }
            for (nm in names(splits)) {
              subg <- splits[[nm]]
              # Only assign fallback when no pairwise is possible in stratum
              if (length(unique(as.character(subg[[factor]]))) <= 1) {
                base <- unique(subg[, c(factor, by_vars), drop = FALSE])
                out <- data.frame(
                  group = as.character(base[[factor]]),
                  letter = rep('a', nrow(base)),
                  stringsAsFactors = FALSE
                )
                if (length(by_vars) > 0) {
                  for (bn in by_vars) {
                    if (bn %in% names(base)) out[[bn]] <- as.character(base[[bn]])
                  }
                }
                assign_list[[length(assign_list) + 1]] <- out
              }
            }
            if (length(assign_list) > 0) {
              cld_df <- do.call(rbind, assign_list)
            } else {
              # Leave as blank to signal computation issue; do not mislead with 'a'
              cld_df <- data.frame(group = character(0), letter = character(0), stringsAsFactors = FALSE)
            }
          } else {
            cld_df <- data.frame(group = character(0), letter = character(0), stringsAsFactors = FALSE)
          }
        }
        if (!is.null(cld_df) && nrow(cld_df) > 0) {
          cld_df$group <- normalizeGroupLabels(cld_df$group, factor)
        }
        tukey_results[[factor]] <- list(test = pw, results = pw_df, letters = cld_df)
      } else {
        tukey_test <- TukeyHSD(anova_model, which = factor)
        tukey_df <- as.data.frame(tukey_test[[factor]])
        tukey_df$comparison <- rownames(tukey_df)
        rownames(tukey_df) <- NULL
        # If too many comparisons, filter to simple differences
        if (shouldFilterComparisons(nrow(tukey_df))) {
          tukey_df <- filterSimpleComparisons(tukey_df)
        }
        if (nrow(tukey_df) > 0) {
          lhs_rhs <- strsplit(tukey_df$comparison, '-')
          lhs <- vapply(lhs_rhs, function(p) trimws(p[[1]]), character(1))
          rhs <- vapply(lhs_rhs, function(p) trimws(p[[2]]), character(1))
          if (identical(factor, 'dbscan_cluster')) {
            lhs <- paste('Cluster', lhs)
            rhs <- paste('Cluster', rhs)
          }
          tukey_df$comparison_pretty <- paste0(factor, ': ', lhs, ' vs ', rhs)
          tukey_df$comparison <- tukey_df$comparison_pretty
        }
        tukey_letters <- generateLabelDf(tukey_test, factor)
        if (!is.null(tukey_letters) && nrow(tukey_letters) > 0) {
          tukey_letters$group <- normalizeGroupLabels(tukey_letters$group, factor)
        }
        tukey_results[[factor]] <- list(test = tukey_test, results = tukey_df, letters = tukey_letters)
      }
    }
    # Ensure the list preserves the exact order of tukey_factors
    tukey_results <- tukey_results[tukey_factors]
    list(results = tukey_results, success = TRUE)
  }, error = function(e) {
    list(error = e$message, success = FALSE)
  })
}

#' Generate label data frame for Tukey results
#' @param TUKEY Tukey test results
#' @param variable Variable name
#' @return Data frame with labels
generateLabelDf <- function(TUKEY, variable) {
  Tukey.levels <- TUKEY[[variable]][, 4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment), ]
  grp <- Tukey.labels$treatment
  if (identical(variable, 'dbscan_cluster')) {
    grp <- paste('Cluster', grp)
  }
  data.frame(letter = Tukey.labels$Letters, group = grp, stringsAsFactors = FALSE)
}


