test_that("PC explain + coach status advance through pipeline", {
  skip_if_not(file.exists(here::here("shiny", "common", "pipeline_coach.R")))
  source(here::here("shiny", "common", "pipeline_coach.R"), local = TRUE)

  st0 <- .sf_coach_status(NULL, NULL, NULL)
  expect_equal(st0$step, 1L)

  ps <- data.frame(unit = 1:3, cultivar = c("A", "B", "C"), mean = 1:3)
  st1 <- .sf_coach_status(ps, NULL, NULL)
  expect_equal(st1$step, 2L)

  loadings <- data.frame(
    trait = c("mean", "auc", "peak"),
    PC1 = c(0.58, 0.58, 0.58),
    PC2 = c(0.01, 0.70, -0.71),
    stringsAsFactors = FALSE
  )
  expect_match(.sf_explain_pc(loadings, "PC1"), "overall curve size")
  expect_match(.sf_explain_pc(loadings, "PC2"), "peak vs integrated AUC|curve shape")

  dr <- list(
    meta = list(success = TRUE),
    loadings = loadings,
    variance = data.frame(component = "PC1", prop_var = 0.99, stringsAsFactors = FALSE)
  )
  st2 <- .sf_coach_status(ps, dr, NULL)
  expect_equal(st2$step, 3L)

  dun <- data.frame(
    contrast = c("B - A", "C - A"),
    p.value = c(0.01, 0.4),
    stringsAsFactors = FALSE
  )
  bl <- list(
    success = TRUE,
    genotype_table = data.frame(level = c("B", "A"), blup = c(1, 0), se = 0.1, rank = 1:2),
    meta = list(trait = "PC1", checks = "A", dunnett = dun)
  )
  st3 <- .sf_coach_status(ps, dr, bl)
  expect_equal(st3$step, 4L)
  expect_equal(nrow(st3$sig), 1L)
})
