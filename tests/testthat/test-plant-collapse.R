# Task 2: ensureTxyColumns + collapsePlantsBySelectedTimes

sf_source_common("txy_parse.R")
sf_source_common("plant_collapse.R")

test_that("ensureTxyColumns parses T:X:Y when geometry cols are absent", {
  raw <- data.frame(
    `T:X:Y` = c("1:2:3", "4:5:6"),
    unit = c("1:2:3", "4:5:6"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- ensureTxyColumns(raw)
  expect_true(all(c("table_no", "pos_x", "pos_y") %in% names(out)))
  expect_equal(out$table_no, c(1L, 4L))
  expect_equal(out$pos_x, c(2L, 5L))
  expect_equal(out$pos_y, c(3L, 6L))
})

test_that("ensureTxyColumns keeps existing table_no/pos_x/pos_y", {
  raw <- data.frame(
    `T:X:Y` = "9:9:9",
    table_no = 1L,
    pos_x = 2L,
    pos_y = 3L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- ensureTxyColumns(raw)
  expect_equal(out$table_no, 1L)
  expect_equal(out$pos_x, 2L)
  expect_equal(out$pos_y, 3L)
})

test_that("ensureTxyColumns parses unit table:x:y when T:X:Y absent (soy-like)", {
  # Soy merged: unit holds geometry; table_no present; pos_x/pos_y absent.
  raw <- data.frame(
    unit = c("51:1:1", "51:1:2", "52:3:4"),
    table_no = c(51L, 51L, 52L),
    cultivar = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  out <- ensureTxyColumns(raw)
  expect_true(all(c("table_no", "pos_x", "pos_y") %in% names(out)))
  expect_equal(out$table_no, c(51L, 51L, 52L))
  expect_equal(out$pos_x, c(1L, 1L, 3L))
  expect_equal(out$pos_y, c(1L, 2L, 4L))
})

test_that("ensureTxyColumns soft-skips when geometry and TXY source absent", {
  raw <- data.frame(
    unit = c("P1", "P2"),
    height = c(1, 2),
    stringsAsFactors = FALSE
  )
  out <- ensureTxyColumns(raw)
  expect_false("pos_x" %in% names(out))
  expect_equal(out$unit, raw$unit)
})

test_that("collapse soft-tolerates missing geometry (non-spatial path)", {
  raw <- data.frame(
    unit = c("P1", "P1", "P2", "P2"),
    cultivar = c("A", "A", "B", "B"),
    treatment = c("T0", "T0", "T0", "T0"),
    height = c(10, 12, 20, 22),
    dbscan_cluster = c(1L, 2L, 1L, 2L),
    timestamp_group = as.POSIXct(
      c("2020-01-01", "2020-01-02", "2020-01-01", "2020-01-02"),
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )
  out <- collapsePlantsBySelectedTimes(
    data = raw,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = "mean"
  )
  expect_equal(nrow(out), 2L)
  expect_true(all(c("unit", "mean") %in% names(out)))
  expect_false("pos_x" %in% names(out))
})

test_that("collapse over two clusters yields one row per plant", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  plants <- dplyr::n_distinct(fix$unit)
  out <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
  expect_equal(nrow(out), plants)
  expect_equal(dplyr::n_distinct(out$unit), plants)
  expect_true(all(c("mean", "auc", "peak") %in% names(out)))
  expect_true(all(c("unit", "V.T.R", "cultivar", "treatment", "table_no", "pos_x", "pos_y") %in% names(out)))
})

test_that("full times vs selected subset differ in collapse metrics", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  full <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = sort(unique(fix$dbscan_cluster)),
    metrics = c("mean", "auc", "peak")
  )
  subset <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
  expect_equal(nrow(full), nrow(subset))
  expect_false(isTRUE(all.equal(full$mean, subset$mean)))
  expect_false(isTRUE(all.equal(full$auc, subset$auc)))
  expect_false(isTRUE(all.equal(full$peak, subset$peak)))
})

test_that("collapsed mini matches committed expected RDS", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  actual <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
  expected_path <- sf_fixture_path("plant_collapsed_mini.expected.rds")
  expect_true(file.exists(expected_path))
  expected <- readRDS(expected_path)
  expect_equal(actual, expected)
})
