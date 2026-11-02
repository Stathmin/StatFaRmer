# Static + smoke suite for soy-325 (no treatment column).
# Skips when local processed RDS is absent (gitignored).

test_that("soy-325 processed manifests exist or skip", {
  root <- sf_repo_root()
  proj <- "project_2025-07-14_b51-91_soy"
  merged <- file.path(root, "data", proj, paste0(proj, "_merged_table.rds"))
  skip_if_not(file.exists(merged), "soy-325 merged RDS not present locally")
  expect_true(file.exists(merged))
})

test_that("soy-325 UI/facet/letters assumptions hold", {
  root <- sf_repo_root()
  proj <- "project_2025-07-14_b51-91_soy"
  merged_path <- file.path(root, "data", proj, paste0(proj, "_merged_table.rds"))
  skip_if_not(file.exists(merged_path), "soy-325 merged RDS not present locally")

  sf_source_common("utils.R")
  sf_source_common("plotting.R")
  sf_source_common("letters_join.R")

  merged <- readRDS(merged_path)
  expect_false("treatment" %in% names(merged) && length(unique(stats::na.omit(merged$treatment))) > 1)

  expect_true(isProjectProcessed(proj))
  sels <- prepareUISelections(merged)
  expect_true(length(sels$anova_factors) > 0)

  # Default facet for 0/1 treatment levels must only use present columns
  n_treat <- length(sels$treatments)
  default_facet <- if (n_treat <= 1L) {
    if ("dbscan_cluster" %in% sels$anova_factors) "~ dbscan_cluster" else "~ ."
  } else {
    "treatment ~ dbscan_cluster"
  }
  fvars <- setdiff(all.vars(as.formula(default_facet)), ".")
  expect_true(all(fvars %in% names(merged)))

  # Facet sanitizer drops missing treatment
  skip_if_not(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 missing")
  p <- ggplot2::ggplot(
    merged[seq_len(min(200L, nrow(merged))), , drop = FALSE],
    ggplot2::aes(x = dbscan_cluster, y = !!rlang::sym(sels$out_variables[1]))
  ) + ggplot2::geom_boxplot()
  p2 <- applyStatfarmerFacet(p, as.formula("treatment ~ dbscan_cluster"), data = merged)
  expect_silent(ggplot2::ggplot_build(p2))

  # Underscore factor letters join
  lev <- as.character(unique(merged$dbscan_cluster))
  letters_df <- data.frame(group = lev, letter = letters[seq_along(lev)], stringsAsFactors = FALSE)
  jl <- joinLettersToData(merged[1:30, , drop = FALSE], letters_df, "dbscan_cluster")
  expect_true("_dbscan_cluster_letter" %in% names(jl$df))
  expect_gt(sum(!is.na(jl$df[["_dbscan_cluster_letter"]])), 0)

  created_temp_cols <- "_dbscan_cluster_letter"
  letter_cols <- created_temp_cols[created_temp_cols %in% names(jl$df)]
  expect_equal(letter_cols, "_dbscan_cluster_letter")
})

test_that("growth grouping prefers columns that exist", {
  root <- sf_repo_root()
  proj <- "project_2025-07-14_b51-91_soy"
  merged_path <- file.path(root, "data", proj, paste0(proj, "_merged_table.rds"))
  skip_if_not(file.exists(merged_path), "soy-325 merged RDS not present locally")
  sf_source_common("utils.R")
  merged <- readRDS(merged_path)
  sels <- prepareUISelections(merged)
  choices <- setdiff(sels$anova_factors, "dbscan_cluster")
  preferred <- intersect(c("treatment", "cultivar", "genotype"), choices)
  selected <- if (length(preferred) > 0L) preferred else head(choices, 1L)
  expect_true(length(selected) >= 1)
  expect_true(all(selected %in% names(merged)))
  expect_false("treatment" %in% selected && !"treatment" %in% names(merged))
})
