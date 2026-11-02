test_that("HortControl GLI/NDVI bin intervals parse correctly", {
  root <- here::here()
  source(file.path(root, "shiny/common/hortcontrol_labels.R"), local = TRUE)
  source(file.path(root, "shiny/common/biology_readout.R"), local = TRUE)
  expect_equal(.sf_spectral_bin_interval("gli_0_0_15_percent"), "[0, 0.15]")
  expect_equal(.sf_spectral_bin_interval("gli_0_2_0_percent"), "[-0.2, 0]")
  expect_equal(.sf_spectral_bin_interval("gli_1_0_2_percent"), "[-1.0, -0.2]")
  expect_equal(.sf_spectral_bin_interval("ndvi_1_0_percent"), "[-1, 0]")
  expect_equal(.sf_hortcontrol_trait_label("gli_0_0_15_percent__auc"), "GLI bin [0, 0.15] · auc")
  expect_equal(.sf_hortcontrol_trait_label("digital_biomass_mm3__last"), "Digital Biomass · last")
  expect_equal(.sf_biology_trait_label("ndvi_average"), "NDVI average")
  # High PC loading on low bins + negative gli_average → NOT greener
  vals <- c(
    gli_0_2_0_percent__auc = 0.28,
    gli_1_0_2_percent__auc = 0.22,
    gli_average__auc = -0.27,
    ndvi_average__auc = -0.24
  )
  expect_false(.sf_biology_high_is_greener(vals))
})

test_that("biology readout builds axes, overlaps, and markdown", {
  root <- here::here()
  source(file.path(root, "shiny/common/pipeline_coach.R"), local = TRUE)
  source(file.path(root, "shiny/common/biology_readout.R"), local = TRUE)

  loadings <- data.frame(
    trait = c("digital_biomass_mm3", "ndvi_average", "saturation_average_percent", "convex_hull_area_mm2"),
    PC1 = c(0.1, -0.2, -0.5, 0.05),
    PC2 = c(0.05, 0.55, 0.1, 0.02),
    PC3 = c(0.6, 0.05, 0.0, 0.55),
    stringsAsFactors = FALSE
  )
  variance <- data.frame(
    component = c("PC1", "PC2", "PC3"),
    prop_var = c(0.3, 0.25, 0.2),
    stringsAsFactors = FALSE
  )
  mk_blup <- function(levels, scores) {
    rk <- rank(-scores, ties.method = "min")
    list(
      success = TRUE,
      genotype_table = data.frame(
        level = levels, blup = scores, se = 0.1, rank = as.integer(rk),
        stringsAsFactors = FALSE
      ),
      meta = list(dunnett = NULL, checks = "G05", trait = "PC")
    )
  }
  set.seed(1)
  lev <- sprintf("G%02d", 1:40)
  blup_by_pc <- list(
    PC1 = mk_blup(lev, rnorm(40)),
    PC2 = mk_blup(lev, rnorm(40)),
    PC3 = mk_blup(lev, rnorm(40))
  )

  r <- buildBiologyReadout(
    loadings, variance, blup_by_pc,
    check = "G05", top_k = 5L, overlap_k = 10L
  )
  expect_true(is.list(r$axes))
  expect_true(nzchar(r$markdown))
  expect_gte(length(r$bullets), 3L)
  expect_match(r$markdown, "Biological readout")
  expect_true(!is.null(r$size_pc))
  expect_true(!is.null(r$green_pc))
})
