test_that("sidebar plant filter + BLUP table filters honor cultivars/treatments", {
  skip_if_not(file.exists(here::here("shiny", "common", "plant_collapse.R")))
  source(here::here("shiny", "common", "plant_collapse.R"), local = TRUE)

  ps <- data.frame(
    unit = 1:6,
    cultivar = c("A", "A", "B", "B", "C", "C"),
    treatment = c("T1", "T2", "T1", "T2", "T1", "T2"),
    mean = rnorm(6),
    stringsAsFactors = FALSE
  )
  ci <- list(cultivars = c("A", "B"), treatments = "T1")
  out <- .sf_filter_plants_by_sidebar(ps, ci)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$cultivar %in% c("A", "B")))
  expect_true(all(out$treatment == "T1"))

  geno <- data.frame(
    level = c("A", "B", "C"),
    blup = c(1, 0, -1),
    se = 0.1,
    rank = 1:3,
    stringsAsFactors = FALSE
  )
  g2 <- .sf_blup_filter_geno_table(geno, ci)
  expect_equal(as.character(g2$level), c("A", "B"))

  trt <- data.frame(
    level = c("T1", "T2"),
    estimate = c(1, 2),
    se = 0.1,
    rank = 1:2,
    stringsAsFactors = FALSE
  )
  t2 <- .sf_blup_filter_trt_table(trt, ci)
  expect_equal(as.character(t2$level), "T1")

  dun <- data.frame(
    contrast = c("B - A", "C - A"),
    estimate = c(0.1, -0.2),
    stringsAsFactors = FALSE
  )
  d2 <- .sf_blup_filter_dunnett(dun, ci)
  expect_true(all(grepl("A|B", d2$contrast)))
  expect_false(any(grepl("^C ", d2$contrast) & !grepl("A|B", d2$contrast)))
})
