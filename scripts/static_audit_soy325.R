#!/usr/bin/env Rscript
# Static audit: soy-325 (no treatment) against Main App / BLUP-DR assumptions.
# Usage (repo root, renv active):
#   Rscript scripts/static_audit_soy325.R

root <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
setwd(root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(rlang)
})

fail <- character(0)
warn <- character(0)
ok <- function(msg) cat("OK  ", msg, "\n", sep = "")
bad <- function(msg) {
  fail <<- c(fail, msg)
  cat("FAIL", msg, "\n", sep = " ")
}
note <- function(msg) {
  warn <<- c(warn, msg)
  cat("WARN", msg, "\n", sep = " ")
}

project_id <- "project_2025-07-14_b51-91_soy"
merged_path <- here("data", project_id, paste0(project_id, "_merged_table.rds"))
groups_path <- here("data", project_id, paste0(project_id, "_vector_of_groups.rds"))

# --- Data presence ---
if (!file.exists(merged_path) || !file.exists(groups_path)) {
  bad(paste("soy-325 processed RDS missing:", merged_path))
  cat("\nSTATIC AUDIT FAILED\n")
  quit(status = 1)
}
merged <- readRDS(merged_path)
ok(sprintf("loaded soy-325: %d rows x %d cols", nrow(merged), ncol(merged)))

has_treatment <- "treatment" %in% names(merged)
n_treat <- if (has_treatment) length(unique(stats::na.omit(merged$treatment))) else 0L
cat(sprintf("INFO treatment column: %s; n_levels=%d\n", has_treatment, n_treat))
if (has_treatment && n_treat > 1L) {
  note("soy-325 unexpectedly has multi-level treatment — audit assumptions may not match")
}

# --- Source helpers ---
source(here("shiny", "common", "utils.R"), local = TRUE)
source(here("shiny", "common", "plotting.R"), local = TRUE)
source(here("shiny", "common", "txy_parse.R"), local = TRUE)
source(here("shiny", "common", "plant_collapse.R"), local = TRUE)
source(here("shiny", "common", "spatial_residualize.R"), local = TRUE)
source(here("shiny", "common", "blup_ranking.R"), local = TRUE)
source(here("shiny", "common", "trait_reduction.R"), local = TRUE)
source(here("shiny", "common", "letters_join.R"), local = TRUE)

# --- Readiness helpers ---
if (!isTRUE(isProjectProcessed(project_id))) {
  bad("isProjectProcessed(soy-325) is FALSE")
} else {
  ok("isProjectProcessed(soy-325)")
}

sels <- prepareUISelections(merged)
if ("treatment" %in% sels$anova_factors && n_treat <= 1L && !has_treatment) {
  bad("prepareUISelections lists treatment without column")
}
if (!has_treatment && length(sels$treatments) > 0L) {
  bad("treatments non-empty without treatment column")
} else {
  ok(sprintf("UI selections: %d anova factors, %d treatments, %d cultivars",
             length(sels$anova_factors), length(sels$treatments), length(sels$cultivars)))
}

# --- Default facet must not require missing treatment ---
n_treat_ui <- length(sels$treatments)
default_facet <- if (n_treat_ui <= 1L) {
  if ("dbscan_cluster" %in% sels$anova_factors) "~ dbscan_cluster" else "~ ."
} else {
  "treatment ~ dbscan_cluster"
}
fml <- as.formula(default_facet)
fvars <- setdiff(all.vars(fml), ".")
missing <- setdiff(fvars, names(merged))
if (length(missing) > 0L) {
  bad(paste("default facet vars missing from data:", paste(missing, collapse = ",")))
} else {
  ok(paste("default facet OK:", default_facet))
}

# --- applyStatfarmerFacet must not die on treatment ~ dbscan_cluster ---
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(merged, ggplot2::aes(x = dbscan_cluster, y = plant_height_averaged_mm)) +
    ggplot2::geom_boxplot()
  p2 <- tryCatch(
    applyStatfarmerFacet(p, as.formula("treatment ~ dbscan_cluster"), ncol = 2L, data = merged),
    error = function(e) e
  )
  if (inherits(p2, "error")) {
    bad(paste("applyStatfarmerFacet crashed:", conditionMessage(p2)))
  } else {
    built <- tryCatch(ggplot2::ggplot_build(p2), error = function(e) e)
    if (inherits(built, "error")) {
      bad(paste("ggplot_build after facet sanitize failed:", conditionMessage(built)))
    } else {
      ok("applyStatfarmerFacet strips missing treatment without combine_vars abort")
    }
  }
} else {
  note("ggplot2 not available; skipped facet build check")
}

# --- Letters join for dbscan_cluster (underscore in name) ---
letters_df <- data.frame(
  group = as.character(unique(as.character(merged$dbscan_cluster))),
  letter = letters[seq_along(unique(as.character(merged$dbscan_cluster)))],
  stringsAsFactors = FALSE
)
jl <- joinLettersToData(merged[seq_len(min(50L, nrow(merged))), , drop = FALSE], letters_df, "dbscan_cluster")
temp_col <- "_dbscan_cluster_letter"
if (!temp_col %in% names(jl$df)) {
  bad("joinLettersToData did not create _dbscan_cluster_letter")
} else if (all(is.na(jl$df[[temp_col]]))) {
  bad("_dbscan_cluster_letter all NA after join")
} else {
  ok("_dbscan_cluster_letter join works")
}
# Simulate reactive_outputs letter_cols selection (no broken [^_] regex)
created_temp_cols <- temp_col
letter_cols <- created_temp_cols[created_temp_cols %in% names(jl$df)]
if (length(letter_cols) != 1L) {
  bad("letter_cols selection failed for underscore factor name")
} else {
  ok("letter_cols selection uses created_temp_cols (underscore-safe)")
}

# --- Static grep: hardcoded treatment defaults in app modules ---
app_files <- list.files(here("shiny", "app"), pattern = "\\.R$", full.names = TRUE)
patterns <- c(
  "selected\\s*=\\s*c\\(\\s*[\"']treatment[\"']",
  "selected\\s*=\\s*[\"']treatment[\"']",
  "value\\s*=\\s*[\"']treatment\\s*~",
  "selected\\s*=\\s*c\\(\\s*[\"']treatment[\"']\\s*,\\s*[\"']cultivar[\"']\\s*\\)"
)
for (f in app_files) {
  lines <- readLines(f, warn = FALSE)
  for (pat in patterns) {
    hits <- grep(pat, lines, perl = TRUE)
    # Allow placeholder text / comments mentioning treatment
    for (h in hits) {
      line <- lines[h]
      if (grepl("placeholder|helpText|#", line)) next
      # ui.R default_facet branch for multi-treat is OK if guarded
      if (grepl("has_multi_treatment|n_treat", paste(lines[max(1, h - 5):h], collapse = "\n"))) next
      note(sprintf("%s:%d hardcoded treatment default: %s", basename(f), h, trimws(line)))
    }
  }
}

# --- Growth grouping selection logic ---
choices <- setdiff(sels$anova_factors, "dbscan_cluster")
preferred <- intersect(c("treatment", "cultivar", "genotype"), choices)
selected <- if (length(preferred) > 0L) preferred else head(choices, 1L)
if (length(selected) < 1L || any(!selected %in% names(merged))) {
  bad(paste("growth grouping selected invalid:", paste(selected, collapse = ",")))
} else {
  ok(paste("growth grouping defaults to:", paste(selected, collapse = ",")))
}

# --- Collapse + BLUP smoke (one trait, subset of clusters) ---
if (!"plant_height_averaged_mm" %in% names(merged)) {
  note("plant_height_averaged_mm missing; skip collapse/BLUP smoke")
} else {
  time_col <- if ("timestamp_group" %in% names(merged)) "timestamp_group" else "dbscan_cluster"
  times <- sort(unique(merged[[time_col]]))
  times <- head(times, 2L)
  sub <- merged[merged[[time_col]] %in% times, , drop = FALSE]
  collapsed <- tryCatch(
    collapsePlantsBySelectedTimes(
      sub,
      trait = "plant_height_averaged_mm",
      time_col = time_col,
      time_values = times,
      metrics = "mean"
    ),
    error = function(e) e
  )
  if (inherits(collapsed, "error")) {
    bad(paste("collapsePlantsBySelectedTimes:", conditionMessage(collapsed)))
  } else {
    ok(sprintf("collapsePlantsBySelectedTimes → %d plant rows", nrow(collapsed)))
    trait_col <- if ("mean" %in% names(collapsed)) "mean" else "plant_height_averaged_mm"
    gvar <- if ("cultivar" %in% names(collapsed)) "cultivar" else "genotype"
    blup <- tryCatch(
      fitBlupModels(
        collapsed,
        trait = trait_col,
        genotype_var = gvar,
        apply_spatial = FALSE,
        top_k = 5L
      ),
      error = function(e) e
    )
    if (inherits(blup, "error")) {
      bad(paste("fitBlupModels:", conditionMessage(blup)))
    } else if (!is.null(blup$ok) && !isTRUE(blup$ok)) {
      bad(paste("fitBlupModels soft-fail:", blup$reason %||% "unknown"))
    } else if (is.null(blup$genotype_table) || nrow(blup$genotype_table) < 1L) {
      bad("fitBlupModels returned empty genotype_table")
    } else {
      ok(sprintf("fitBlupModels genotype_table n=%d", nrow(blup$genotype_table)))
    }
  }
}

cat("\n==== SUMMARY ====\n")
cat(sprintf("FAILS: %d\n", length(fail)))
cat(sprintf("WARNS: %d\n", length(warn)))
if (length(fail) > 0L) {
  cat(paste(" -", fail), sep = "\n")
  cat("\nSTATIC AUDIT FAILED\n")
  quit(status = 1)
}
cat("\nSTATIC AUDIT PASSED\n")
