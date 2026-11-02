test_that("across-time raw mean posthoc is hard-blocked", {
  root <- here::here()
  source(file.path(root, "shiny/common/trait_followup.R"), local = TRUE)

  expect_error(
    .sf_forbid_across_time_raw_mean_posthoc("unit-test"),
    "antipattern blocked"
  )
  expect_error(
    .sf_plant_trait_means(data.frame(x = 1), "x"),
    "antipattern blocked"
  )
  expect_error(
    .sf_plant_trait_means(),
    "do not average raw phenotype"
  )
})

test_that("single-time slice refuses multi-time averaging", {
  root <- here::here()
  source(file.path(root, "shiny/common/trait_followup.R"), local = TRUE)

  d <- data.frame(
    unit = rep(1:4, each = 3),
    cultivar = rep(c("A", "B"), each = 6),
    dbscan_cluster = rep(1:3, times = 4),
    y = rnorm(12),
    stringsAsFactors = FALSE
  )
  sl <- .sf_slice_single_time(d, time_value = NULL)
  expect_equal(sl$time_value, "3")
  expect_equal(length(unique(sl$data$dbscan_cluster)), 1L)
  expect_match(sl$note, "Not an across-time mean|Single time slice")

  sl2 <- .sf_slice_single_time(d, time_value = "1")
  expect_equal(sl2$time_value, "1")
  expect_true(all(as.character(sl2$data$dbscan_cluster) == "1"))
})

test_that("trait follow-up builds curve/violin and single-time Tukey letters", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("emmeans")
  root <- here::here()
  source(file.path(root, "shiny/common/plotting.R"), local = TRUE)
  source(file.path(root, "shiny/common/letters_join.R"), local = TRUE)
  source(file.path(root, "shiny/common/letters_utils.R"), local = TRUE)
  source(file.path(root, "shiny/common/letters_normalize.R"), local = TRUE)
  source(file.path(root, "shiny/common/tukey_filter.R"), local = TRUE)
  source(file.path(root, "shiny/common/letters_emmeans.R"), local = TRUE)
  source(file.path(root, "shiny/common/letters_pvals.R"), local = TRUE)
  source(file.path(root, "shiny/common/tukey.R"), local = TRUE)
  source(file.path(root, "shiny/common/trait_followup.R"), local = TRUE)

  set.seed(1)
  d <- data.frame(
    unit = rep(1:12, each = 3),
    cultivar = rep(rep(c("A", "B", "C"), each = 4), each = 3),
    dbscan_cluster = rep(1:3, times = 12),
    digital_biomass_mm3 = NA_real_,
    stringsAsFactors = FALSE
  )
  d$digital_biomass_mm3 <- as.numeric(factor(d$cultivar)) * 5 +
    d$dbscan_cluster + rnorm(nrow(d), 0, 0.3)

  # Default = last time slice (cluster 3), not mean across 1:3
  b <- buildTraitFollowupPlots(d, "digital_biomass_mm3", group = "cultivar", run_tukey = TRUE)
  expect_s3_class(b$curve, "ggplot")
  expect_s3_class(b$violin, "ggplot")
  expect_equal(as.character(b$time_value), "3")
  expect_match(b$slice_note, "Single time|Not an across-time")
  expect_true(isTRUE(b$tukey$success))
  expect_true(nrow(b$tukey$letters) >= 2L)
  expect_equal(
    length(unique(as.character(b$tukey$slice$data$dbscan_cluster))),
    1L
  )
  # plant_df rows = plants at that time only (≤ n units)
  expect_lte(nrow(b$tukey$plant_df), length(unique(d$unit)))

  b1 <- buildTraitFollowupPlots(
    d, "digital_biomass_mm3",
    group = "cultivar", time_value = "1", run_tukey = TRUE
  )
  expect_equal(as.character(b1$time_value), "1")
  expect_true(all(as.character(b1$tukey$slice$data$dbscan_cluster) == "1"))
})
