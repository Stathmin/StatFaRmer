test_that("collapse supports growth_rate and last", {
  root <- here::here()
  source(file.path(root, "shiny/common/txy_parse.R"), local = TRUE)
  source(file.path(root, "shiny/common/plant_collapse.R"), local = TRUE)

  raw <- data.frame(
    unit = c("P1", "P1", "P1", "P2", "P2", "P2"),
    cultivar = c("A", "A", "A", "B", "B", "B"),
    treatment = "T0",
    height = c(10, 20, 40, 5, 8, 9),
    dbscan_cluster = c(1L, 2L, 3L, 1L, 2L, 3L),
    timestamp_group = as.POSIXct(
      c("2020-01-01", "2020-01-02", "2020-01-04",
        "2020-01-01", "2020-01-02", "2020-01-04"),
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )
  out <- collapsePlantsBySelectedTimes(
    raw, "height", "dbscan_cluster", c(1L, 2L, 3L),
    metrics = c("auc", "peak", "growth_rate", "last")
  )
  expect_equal(nrow(out), 2L)
  expect_true(all(c("auc", "peak", "growth_rate", "last") %in% names(out)))
  expect_equal(out$last[out$unit == "P1"], 40)
  expect_equal(out$peak[out$unit == "P1"], 40)
  # (40-10) / 3 days
  expect_equal(out$growth_rate[out$unit == "P1"], 10)
})

test_that("collapseManyTraits yields trait__metric columns and assert passes", {
  root <- here::here()
  source(file.path(root, "shiny/common/txy_parse.R"), local = TRUE)
  source(file.path(root, "shiny/common/plant_collapse.R"), local = TRUE)

  raw <- data.frame(
    unit = rep(c("P1", "P2"), each = 3),
    cultivar = rep(c("A", "B"), each = 3),
    treatment = "T0",
    digital_biomass_mm3 = c(1, 2, 4, 1, 1.5, 2),
    ndvi_average = c(0.4, 0.5, 0.55, 0.3, 0.35, 0.4),
    dbscan_cluster = rep(1:3, 2),
    timestamp_group = as.POSIXct(
      rep(c("2020-01-01", "2020-01-02", "2020-01-03"), 2),
      tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )
  wide <- collapseManyTraitsBySelectedTimes(
    raw,
    traits = c("digital_biomass_mm3", "ndvi_average"),
    time_col = "dbscan_cluster",
    time_values = 1:3,
    metrics = c("auc", "last")
  )
  expect_true(all(c(
    "digital_biomass_mm3__auc", "digital_biomass_mm3__last",
    "ndvi_average__auc", "ndvi_average__last"
  ) %in% names(wide)))
  mc <- attr(wide, "metric_cols")
  expect_silent(.sf_dr_assert_named_temporal_metrics(mc))
})

test_that("unlabeled raw phenotype columns fail DR temporal assert", {
  root <- here::here()
  source(file.path(root, "shiny/common/plant_collapse.R"), local = TRUE)
  expect_error(
    .sf_dr_assert_named_temporal_metrics(c("digital_biomass_mm3", "ndvi_average")),
    "named temporal metrics"
  )
  expect_error(
    .sf_dr_assert_named_temporal_metrics(c("digital_biomass_mm3__auc", "ndvi_average")),
    "named temporal metrics"
  )
})

test_that("reduceTraitsPCA refuses unlabeled raw traits when required", {
  root <- here::here()
  source(file.path(root, "shiny/common/plant_collapse.R"), local = TRUE)
  source(file.path(root, "shiny/common/trait_reduction.R"), local = TRUE)

  plant <- data.frame(
    unit = 1:10,
    cultivar = rep(LETTERS[1:5], each = 2),
    treatment = "ctrl",
    digital_biomass_mm3 = rnorm(10, 10),
    ndvi_average = rnorm(10, 0.5),
    stringsAsFactors = FALSE
  )
  bad <- reduceTraitsPCA(
    plant,
    trait_cols = c("digital_biomass_mm3", "ndvi_average"),
    require_named_temporal_metrics = TRUE
  )
  expect_false(isTRUE(bad$meta$success))
  expect_match(bad$meta$reason, "Named temporal metrics|unlabeled")

  plant2 <- plant
  names(plant2)[names(plant2) == "digital_biomass_mm3"] <- "digital_biomass_mm3__auc"
  names(plant2)[names(plant2) == "ndvi_average"] <- "ndvi_average__last"
  ok <- reduceTraitsPCA(
    plant2,
    trait_cols = c("digital_biomass_mm3__auc", "ndvi_average__last"),
    require_named_temporal_metrics = TRUE
  )
  expect_true(isTRUE(ok$meta$success))
})
