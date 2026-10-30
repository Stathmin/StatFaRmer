# BLUP / Trait-DR plots must use shared plotting.R helpers (ggplot + theme).

test_that("sidecar plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")

  plotting <- here::here("shiny", "common", "plotting.R")
  skip_if_not(file.exists(plotting))
  source(plotting, local = TRUE)

  g <- data.frame(
    level = paste0("G", 1:5),
    blup = c(1, -0.2, 0.5, 0.1, -0.8),
    se = 0.1,
    rank = 1:5,
    stringsAsFactors = FALSE
  )
  expect_s3_class(createBlupRankPlot(g, checks = "G1"), "ggplot")

  d <- data.frame(
    contrast = c("G2 - G1", "G3 - G1"),
    estimate = c(0.2, -0.3),
    SE = c(0.1, 0.1),
    p.value = c(0.2, 0.01),
    stringsAsFactors = FALSE
  )
  attr(d, "dunnett_ref") <- "G1"
  expect_s3_class(createBlupDunnettPlot(d), "ggplot")
  expect_s3_class(
    createBlupTreatmentPlot(
      data.frame(level = c("T1", "T2"), estimate = c(1, 2), se = c(0.1, 0.2), rank = 1:2)
    ),
    "ggplot"
  )

  v <- data.frame(
    component = c("PC1", "PC2"),
    prop_var = c(0.6, 0.3),
    std_dev = c(1, 0.7),
    cum_prop = c(0.6, 0.9),
    stringsAsFactors = FALSE
  )
  expect_s3_class(createDrScreePlot(v), "ggplot")

  pl <- data.frame(
    unit = 1:10,
    cultivar = rep(c("A", "B"), 5),
    PC1 = rnorm(10),
    PC2 = rnorm(10),
    stringsAsFactors = FALSE
  )
  expect_s3_class(createDrScoresPlot(pl), "ggplot")

  ld <- data.frame(
    trait = c("mean", "auc", "peak"),
    PC1 = c(0.7, 0.2, -0.1),
    PC2 = c(0.1, 0.8, 0.3),
    stringsAsFactors = FALSE
  )
  expect_s3_class(createDrLoadingsPlot(ld), "ggplot")
  expect_s3_class(statfarmerEmptyPlot("idle"), "ggplot")
})
