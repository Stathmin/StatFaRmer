# Assumption checks utilities

library(tidyverse)
library(car)

#' Perform assumption checks for ANOVA/mixed models
#' @param model aov or lmerMod model
#' @param data original data
#' @param response response variable name
#' @param main_factor main factor for variance test (first factor on RHS)
#' @return tibble with shapiro and levene/bartlett results
performAssumptionChecks <- function(model, data, response, main_factor) {
  is_lmer <- inherits(model, 'lmerMod')
  res <- residuals(model)
  shapiro_p <- tryCatch({
    shapiro.test(res)$p.value
  }, error = function(e) NA_real_)
  levene_p <- if (!is_lmer) {
    tryCatch({
      if (!is.null(main_factor) && main_factor %in% names(data)) {
        car::leveneTest(data[[response]] ~ as.factor(data[[main_factor]]), data = data)["Pr(>F)"][1, 1]
      } else NA_real_
    }, error = function(e) NA_real_)
  } else NA_real_
  bartlett_p <- if (!is_lmer) {
    tryCatch({
      if (!is.null(main_factor) && main_factor %in% names(data)) {
        bartlett.test(data[[response]] ~ as.factor(data[[main_factor]]), data = data)$p.value
      } else NA_real_
    }, error = function(e) NA_real_)
  } else NA_real_
  result <- tibble::tibble(
    check = c("Shapiro-Wilk (normality)", "Levene (homogeneity)", "Bartlett (homogeneity)"),
    p_value = c(shapiro_p, levene_p, bartlett_p),
    pass = dplyr::case_when(
      is.na(p_value) ~ NA,
      TRUE ~ p_value > 0.05
    )
  )
  if (is_lmer) {
    result$note <- c(
      "Tests conditional residuals",
      "Not applicable for mixed models",
      "Not applicable for mixed models"
    )
  } else {
    result$note <- rep(NA_character_, 3)
  }
  result
}


