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
    tukey_results <- list()
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
        if (length(other) > 0) {
          spec_str <- paste(factor, '|', paste(other, collapse = ' + '))
          emm <- emmeans::emmeans(anova_model, as.formula(paste('~', spec_str)))
        } else {
          emm <- emmeans::emmeans(anova_model, as.formula(paste('~', factor)))
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
        cld_df <- tryCatch({
          pvals <- setNames(pw_df$p.adj, pw_df$contrast)
          letters_result <- multcompView::multcompLetters(pvals)
          factor_levels <- as.character(unique(summary(emm)[[factor]]))
          letters_names <- names(letters_result$Letters)
          letters_names <- trimws(letters_names)
          stripped_names <- sub(paste0('^', factor, '(.*)$'), '\\1', letters_names)
          names(letters_result$Letters) <- stripped_names
          letters_vec <- letters_result$Letters[factor_levels]
          data.frame(group = factor_levels, letter = as.character(letters_vec), stringsAsFactors = FALSE)
        }, error = function(e) {
          data.frame(group = character(0), letter = character(0))
        })
        tukey_results[[factor]] <- list(test = pw, results = pw_df, letters = cld_df)
      } else {
        tukey_test <- TukeyHSD(anova_model, which = factor)
        tukey_df <- as.data.frame(tukey_test[[factor]])
        tukey_df$comparison <- rownames(tukey_df)
        rownames(tukey_df) <- NULL
        parts <- strsplit(tukey_df$comparison, '-')
        keep_mask <- vapply(parts, function(p) {
          if (length(p) != 2) return(FALSE)
          lhs <- trimws(p[[1]]); rhs <- trimws(p[[2]])
          lhs_levels <- strsplit(lhs, ':')[[1]]
          rhs_levels <- strsplit(rhs, ':')[[1]]
          if (length(lhs_levels) != length(rhs_levels)) return(FALSE)
          sum(lhs_levels != rhs_levels) == 1
        }, logical(1))
        tukey_df <- tukey_df[keep_mask, , drop = FALSE]
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
        tukey_results[[factor]] <- list(test = tukey_test, results = tukey_df, letters = tukey_letters)
      }
    }
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


