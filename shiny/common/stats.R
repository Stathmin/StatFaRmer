# StatFaRmer Statistical Functions
# Iteration 3: Modular Shiny
# Iteration 7: Time-Aware Mixed Models

# Load required libraries
library(tidyverse)
library(broom)
library(emmeans)
library(multcomp)
library(multcompView)
library(moments)
library(flextable)
library(car)
library(lme4)

# Source split modules to keep this file high-level
try({
  if (requireNamespace('here', quietly = TRUE)) {
    src <- function(f) if (file.exists(here::here('shiny','common', f))) source(here::here('shiny','common', f))
    src('anova.R'); src('spline.R'); src('tukey.R'); src('assumptions.R'); src('model_resolver.R')
  }
}, silent = TRUE)

# =============================================================================
# LEGACY FUNCTIONS MOVED TO MODULES
# =============================================================================
# analyzeFactorCardinality -> anova.R
# performANOVA -> model_resolver.R (delegated below)
# performAssumptionChecks -> assumptions.R
# performTukey, generateLabelDf -> tukey.R
# formatANOVATable, formatDescTable, formatTukeyTable, formatGroupLettersTable -> tables.R
# createANOVAPlot -> plotting.R

# -----------------------------------------------------------------------------
# Override performANOVA to delegate to resolver (kept at end to mask older def)
# -----------------------------------------------------------------------------
try({
performANOVA <- function(data, formula, out_variable, force_method = NULL) {
    if (exists('resolveModel', mode = 'function')) {
      return(resolveModel(data, formula, out_variable, force_method))
    }
    # fallback: basic aov
    res <- try(aov(as.formula(formula), data = data), silent = TRUE)
    if (inherits(res, 'try-error')) return(list(error = as.character(res), success = FALSE))
    list(
      model = res,
      model_type = 'aov',
      anova_table = broom::tidy(res),
      desc_stats = tibble::tibble(),
      cardinality = integer(0),
      blocking_factors = character(0),
      success = TRUE
    )
  }
}, silent = TRUE)
