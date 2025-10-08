# Table formatting utilities

library(flextable)
library(tidyverse)

#' Format ANOVA table for display
formatANOVATable <- function(anova_table) {
  if (is.null(anova_table) || nrow(anova_table) == 0) {
    return(flextable(data.frame(Message = "No ANOVA results available")))
  }
  anova_table$significance <- dplyr::case_when(
    anova_table$p.value < 0.001 ~ '***',
    anova_table$p.value < 0.01 ~ '**',
    anova_table$p.value < 0.05 ~ '*',
    anova_table$p.value < 0.1 ~ '.',
    TRUE ~ ' '
  )
  anova_table$p.value <- format(anova_table$p.value, digits = 4, scientific = TRUE)
  flextable(anova_table) %>%
    set_header_labels(
      term = "Term",
      df = "DF",
      sumsq = "Sum Sq",
      meansq = "Mean Sq",
      statistic = "F value",
      p.value = "Pr(>F)",
      significance = "Signif."
    ) %>%
    theme_zebra() %>%
    autofit()
}

#' Format descriptive statistics table
formatDescTable <- function(desc_stats) {
  if (is.null(desc_stats) || nrow(desc_stats) == 0) {
    return(flextable(data.frame(Message = "No descriptive statistics available")))
  }
  numeric_cols <- sapply(desc_stats, is.numeric)
  desc_stats[numeric_cols] <- round(desc_stats[numeric_cols], 3)
  flextable(desc_stats) %>% theme_zebra() %>% autofit()
}

#' Format Tukey results table
formatTukeyTable <- function(tukey_results) {
  if (is.null(tukey_results) || length(tukey_results) == 0) {
    return(flextable(data.frame(Message = "No Tukey results available")))
  }
  all_results <- data.frame()
  for (factor_name in names(tukey_results)) {
    factor_results <- tukey_results[[factor_name]]$results
    factor_results$factor <- factor_name
    all_results <- rbind(all_results, factor_results)
  }
  if (nrow(all_results) == 0) {
    return(flextable(data.frame(Message = "No Tukey results available")))
  }
  all_results$significance <- dplyr::case_when(
    all_results$p.adj < 0.001 ~ '***',
    all_results$p.adj < 0.01 ~ '**',
    all_results$p.adj < 0.05 ~ '*',
    all_results$p.adj < 0.1 ~ '.',
    TRUE ~ ' '
  )
  all_results$p.adj <- format(all_results$p.adj, digits = 4, scientific = TRUE)
  if ('comparison_pretty' %in% names(all_results)) {
    all_results$comparison <- all_results$comparison_pretty
    all_results$comparison_pretty <- NULL
  }
  flextable(all_results) %>%
    set_header_labels(
      factor = "Factor",
      comparison = "Comparison",
      diff = "Difference",
      lwr = "Lower CI",
      upr = "Upper CI",
      p.adj = "Adj. P-value",
      significance = "Signif."
    ) %>%
    theme_zebra() %>%
    autofit()
}

#' Format group letters table
formatGroupLettersTable <- function(tukey_results) {
  if (is.null(tukey_results) || length(tukey_results) == 0) {
    return(flextable(data.frame(Message = "No group letters available")))
  }
  all_letters <- data.frame()
  for (factor_name in names(tukey_results)) {
    factor_letters <- tukey_results[[factor_name]]$letters
    factor_letters$factor <- factor_name
    all_letters <- rbind(all_letters, factor_letters)
  }
  if (nrow(all_letters) == 0) {
    return(flextable(data.frame(Message = "No group letters available")))
  }
  flextable(all_letters) %>%
    set_header_labels(
      factor = "Factor",
      group = "Group",
      letter = "Letter"
    ) %>%
    theme_zebra() %>%
    autofit()
}


