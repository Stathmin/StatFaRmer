# Task 4: residualizeSpatial (per-table mgcv smooths)

sf_source_common("txy_parse.R")
sf_source_common("spatial_residualize.R")

test_that("with X/Y, spatial residual column exists and variance drops", {
  skip_if_not_installed("mgcv")
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  out <- residualizeSpatial(fix, y_col = "height")
  expect_type(out, "list")
  expect_true(all(c("data", "meta") %in% names(out)))
  expect_true(isTRUE(out$meta$spatial_applied))
  resid_col <- "height_spat_resid"
  expect_true(resid_col %in% names(out$data))
  expect_equal(length(out$data[[resid_col]]), nrow(fix))
  expect_lt(stats::var(out$data[[resid_col]], na.rm = TRUE), stats::var(fix$height, na.rm = TRUE))
  expect_type(out$meta$formula, "character")
  expect_true(isTRUE(out$meta$converged) || isFALSE(out$meta$converged))
  expect_true(length(out$meta$per_table_n) >= 1L)
})

test_that("without X/Y, passthrough identity residuals", {
  raw <- data.frame(
    height = c(1, 2, 3, 4),
    table_no = c(1L, 1L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  out <- residualizeSpatial(raw, y_col = "height")
  expect_false(isTRUE(out$meta$spatial_applied))
  expect_equal(out$data$height_spat_resid, raw$height)
})

test_that("spatial meta mini matches committed expected RDS", {
  skip_if_not_installed("mgcv")
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  actual <- residualizeSpatial(fix, y_col = "height")$meta
  expected_path <- sf_fixture_path("spatial_meta_mini.expected.rds")
  expect_true(file.exists(expected_path))
  expected <- readRDS(expected_path)
  expect_equal(actual$spatial_applied, expected$spatial_applied)
  expect_equal(actual$per_table_n, expected$per_table_n)
  expect_equal(actual$formula, expected$formula)
  expect_equal(actual$converged, expected$converged)
  expect_equal(actual$resid_col, expected$resid_col)
  expect_equal(actual$k_x, expected$k_x)
  expect_equal(actual$k_y, expected$k_y)
})
