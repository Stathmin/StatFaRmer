test_that("plant_curve_mini fixture has plant-time structure", {
  fix <- readRDS(testthat::test_path("fixtures", "plant_curve_mini.rds"))
  expect_true(all(c("unit", "cultivar") %in% names(fix)))
  expect_gte(dplyr::n_distinct(fix$unit), 8)
  expect_gte(dplyr::n_distinct(fix$dbscan_cluster), 2)
})
