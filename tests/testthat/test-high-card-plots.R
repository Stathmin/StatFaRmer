test_that("high-cardinality plots avoid 325-way legends/violins/Tukey", {
  skip_if_not_installed("ggplot2")
  root <- here::here()
  source(file.path(root, "shiny/common/plotting.R"), local = TRUE)
  source(file.path(root, "shiny/common/trait_followup.R"), local = TRUE)

  set.seed(2)
  n_g <- 60L
  d <- data.frame(
    unit = rep(seq_len(n_g), each = 3),
    cultivar = rep(sprintf("G%02d", seq_len(n_g)), each = 3),
    dbscan_cluster = rep(1:3, times = n_g),
    y = rnorm(n_g * 3, mean = rep(seq_len(n_g), each = 3) / 10),
    stringsAsFactors = FALSE
  )

  curve <- createTraitGrowthCurvePlot(d, "y", group = "cultivar")
  expect_s3_class(curve, "ggplot")
  expect_true(grepl("ribbon|genotypes|sample", curve$labels$subtitle, ignore.case = TRUE))

  strip <- createTraitViolinPlot(d[d$dbscan_cluster == 3, ], "y", group = "cultivar")
  expect_s3_class(strip, "ggplot")
  expect_match(paste(strip$labels$subtitle, collapse = " "), "Top|bottom|genotypes")

  rk <- data.frame(
    level = sprintf("G%02d", seq_len(n_g)),
    blup = rnorm(n_g),
    se = 0.1,
    rank = seq_len(n_g),
    stringsAsFactors = FALSE
  )
  p_rank <- createBlupRankPlot(rk, checks = "G30", keep_ends = 10L)
  expect_lte(attr(p_rank, "sf_n_rows_plotted"), 21L)
  expect_gte(attr(p_rank, "sf_n_rows_plotted"), 20L)

  tuk <- computeTraitTukeyLetters(d, "y", group = "cultivar", time_value = "3")
  expect_false(isTRUE(tuk$success))
  expect_match(tuk$reason, "Tukey CLD skipped")
})
