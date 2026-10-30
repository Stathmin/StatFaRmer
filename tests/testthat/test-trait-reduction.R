# Task 7: plant × traits PCA (reduceTraitsPCA)

sf_source_common("trait_reduction.R")

.sf_wide_plants <- function(n = 12L) {
  set.seed(7)
  data.frame(
    unit = paste0("P", seq_len(n)),
    cultivar = rep(paste0("G", 1:4), length.out = n),
    treatment = rep(c("T0", "T1"), length.out = n),
    `repeat` = rep(1:2, length.out = n),
    table_no = rep(1:3, length.out = n),
    height = rnorm(n, 20, 3),
    area = rnorm(n, 100, 15),
    volume = rnorm(n, 50, 8),
    ndvi = rnorm(n, 0.6, 0.1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("PCA scores have one row per plant", {
  plant_df <- .sf_wide_plants()
  traits <- c("height", "area", "volume", "ndvi")
  out <- reduceTraitsPCA(plant_df, trait_cols = traits, n_comp = 2L)
  expect_equal(nrow(out$plant_df_with_pcs), nrow(plant_df))
  pc_cols <- grep("^PC[0-9]+$", names(out$plant_df_with_pcs), value = TRUE)
  expect_true(length(pc_cols) >= 1L)
  scores <- out$plant_df_with_pcs[, pc_cols, drop = FALSE]
  expect_equal(nrow(scores), nrow(plant_df))
})

test_that("cultivar and design factors are never PCA variables", {
  plant_df <- .sf_wide_plants()
  # Even if a caller mistakenly passes cultivar, it must be excluded
  out <- reduceTraitsPCA(
    plant_df,
    trait_cols = c("height", "area", "volume", "cultivar"),
    n_comp = 2L
  )
  expect_false("cultivar" %in% out$meta$trait_cols)
  expect_false("treatment" %in% out$meta$trait_cols)
  expect_true(all(out$meta$trait_cols %in% c("height", "area", "volume", "ndvi")))
  expect_false("cultivar" %in% rownames(as.matrix(out$loadings)))
})

test_that("PC columns bind onto the same plant rows", {
  plant_df <- .sf_wide_plants()
  out <- reduceTraitsPCA(
    plant_df,
    trait_cols = c("height", "area", "volume", "ndvi"),
    n_comp = 3L
  )
  pcs <- out$plant_df_with_pcs
  expect_true(all(c("unit", "cultivar", "PC1") %in% names(pcs)))
  expect_equal(pcs$unit, plant_df$unit)
  expect_equal(pcs$cultivar, plant_df$cultivar)
  expect_true(is.numeric(pcs$PC1))
  expect_true(is.numeric(pcs$PC2))
})

test_that("default trait_cols exclude design IDs", {
  plant_df <- .sf_wide_plants()
  out <- reduceTraitsPCA(plant_df, trait_cols = NULL, n_comp = 2L)
  expect_true(all(c("height", "area", "volume", "ndvi") %in% out$meta$trait_cols))
  expect_false(any(c("cultivar", "treatment", "table_no", "unit") %in%
    out$meta$trait_cols))
})

test_that("stale PC frame is not preferred after plant summaries change", {
  plant_df <- .sf_wide_plants()
  out <- reduceTraitsPCA(
    plant_df,
    trait_cols = c("height", "area", "volume", "ndvi"),
    n_comp = 2L
  )
  pcs <- out$plant_df_with_pcs
  expect_true(.sf_dr_pcs_aligned(plant_df, pcs))
  expect_identical(.sf_dr_plants_for_blup(plant_df, pcs), pcs)

  refreshed <- plant_df
  refreshed$height <- refreshed$height + 1
  expect_false(.sf_dr_pcs_aligned(refreshed, pcs))
  expect_identical(.sf_dr_plants_for_blup(refreshed, pcs), refreshed)

  fewer <- plant_df[-1, , drop = FALSE]
  expect_false(.sf_dr_pcs_aligned(fewer, pcs))
  expect_identical(.sf_dr_plants_for_blup(fewer, pcs), fewer)
})

