# Task 3: Growth Summaries curve bridge — selected times contract
# Growth Summaries must call collapsePlantsBySelectedTimes with the user's
# selected dbscan_cluster values, not the full unfiltered time series.

sf_source_common("txy_parse.R")
sf_source_common("plant_collapse.R")

test_that("Growth Summaries path: selected clusters collapse ≠ all clusters", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  all_clusters <- sort(unique(fix$dbscan_cluster))
  selected <- all_clusters[seq_len(min(2L, length(all_clusters)))]
  expect_lt(length(selected), length(all_clusters))

  # Contract Growth Summaries must use (same helper + args shape)
  full <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = all_clusters,
    metrics = c("mean", "auc", "peak")
  )
  subset <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = selected,
    metrics = c("mean", "auc", "peak")
  )

  expect_equal(nrow(subset), dplyr::n_distinct(fix$unit))
  expect_equal(nrow(full), nrow(subset))
  expect_false(isTRUE(all.equal(full$mean, subset$mean)))
  expect_false(isTRUE(all.equal(full$auc, subset$auc)))
  expect_false(isTRUE(all.equal(full$peak, subset$peak)))
})

test_that("Growth Summaries path: selected-time collapse matches expected RDS", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  actual <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
  expected <- readRDS(sf_fixture_path("plant_collapsed_mini.expected.rds"))
  expect_equal(actual, expected)
})

test_that("Growth Summaries path: collapsing all times is not the selected-times result", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  selected_only <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = "mean"
  )
  # Simulating the pre-Task-3 bug: pass every cluster as if using full merged_table
  buggy_full_series <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = sort(unique(fix$dbscan_cluster)),
    metrics = "mean"
  )
  expect_false(isTRUE(all.equal(selected_only$mean, buggy_full_series$mean)))
})
