# Task 5: fitBlupModels (genotype BLUP helper)

sf_source_common("txy_parse.R")
sf_source_common("spatial_residualize.R")
sf_source_common("plant_collapse.R")
sf_source_common("blup_ranking.R")

.sf_collapsed_mini <- function() {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
}

test_that("genotype ranks length equals n genotypes on collapsed mini", {
  skip_if_not_installed("lme4")
  collapsed <- .sf_collapsed_mini()
  n_g <- length(unique(collapsed$cultivar))
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE
  )
  expect_true(isTRUE(out$success) || is.null(out$success) || !isFALSE(out$success))
  expect_equal(nrow(out$genotype_table), n_g)
  expect_equal(length(out$genotype_table$rank), n_g)
  expect_equal(sort(out$genotype_table$rank), seq_len(n_g))
  expect_true(all(c("level", "blup", "rank") %in% names(out$genotype_table)))
})

test_that("1-level treatment yields empty treatment_table + message", {
  skip_if_not_installed("lme4")
  collapsed <- .sf_collapsed_mini()
  expect_equal(length(unique(collapsed$treatment)), 1L)
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE
  )
  expect_equal(nrow(out$treatment_table), 0L)
  expect_true(is.character(out$meta$treatment_message) || is.character(out$meta$message))
  msg <- if (!is.null(out$meta$treatment_message)) {
    out$meta$treatment_message
  } else {
    out$meta$message
  }
  expect_match(msg, "treatment|level", ignore.case = TRUE)
})

test_that("multi-time raw input is blocked with success=FALSE", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  out <- fitBlupModels(
    data = fix,
    trait = "height",
    apply_spatial = FALSE
  )
  expect_false(isTRUE(out$success))
  expect_true(is.character(out$reason) || is.character(out$meta$reason))
  expect_true(isTRUE(out$meta$n_time_per_plant > 1) ||
    isTRUE(max(out$meta$n_time_per_plant, na.rm = TRUE) > 1))
})

test_that("genotype_table matches committed expected RDS", {
  skip_if_not_installed("lme4")
  collapsed <- .sf_collapsed_mini()
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE,
    top_k = 20L
  )
  expected_path <- sf_fixture_path("blup_geno_mini.expected.rds")
  expect_true(file.exists(expected_path))
  expected <- readRDS(expected_path)
  actual <- out$genotype_table
  expect_equal(actual$level, expected$level)
  expect_equal(actual$rank, expected$rank)
  expect_equal(actual$blup, expected$blup, tolerance = 1e-6)
  if ("se" %in% names(expected) && "se" %in% names(actual)) {
    expect_equal(actual$se, expected$se, tolerance = 1e-6)
  }
})

test_that("checks non-NULL yields Dunnett contrast data.frame", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  collapsed <- .sf_collapsed_mini()
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE,
    checks = "cvA"
  )
  expect_true(isTRUE(out$success))
  expect_true(is.data.frame(out$meta$dunnett))
  expect_gt(nrow(out$meta$dunnett), 0L)
  expect_true(any(grepl("p.?value|p.value|pvalue", names(out$meta$dunnett), ignore.case = TRUE)))
  expect_null(out$meta$dunnett_message)
})

test_that("skipped spatial keeps original trait in formula and meta$y_col", {
  skip_if_not_installed("lme4")
  # No X/Y → residualizeSpatial identity; must not rename response to *_spat_resid
  raw <- data.frame(
    mean = c(10, 12, 14, 16, 11, 13, 15, 17),
    cultivar = rep(c("cvA", "cvB", "cvC", "cvD"), each = 2L),
    treatment = "ctrl",
    `repeat` = rep(1:2, 4L),
    table_no = rep(1:2, each = 4L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- fitBlupModels(
    data = raw,
    trait = "mean",
    apply_spatial = TRUE
  )
  expect_true(isTRUE(out$success))
  expect_false(isTRUE(out$meta$spatial_applied))
  expect_equal(out$meta$y_col, "mean")
  expect_match(out$formula, "^mean ~")
  expect_false(grepl("_spat_resid", out$formula, fixed = TRUE))
})
