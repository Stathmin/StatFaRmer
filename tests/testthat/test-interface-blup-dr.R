# Task 10: Full static interface gate for BLUP/DR (no Shiny launch).
# Freezes public helper signatures, return shapes, and mini-fixture contracts.

sf_source_common("txy_parse.R")
sf_source_common("spatial_residualize.R")
sf_source_common("plant_collapse.R")
sf_source_common("blup_ranking.R")
sf_source_common("trait_reduction.R")

.sf_iface_collapsed_mini <- function() {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = c(1L, 2L),
    metrics = c("mean", "auc", "peak")
  )
}

.sf_iface_wide_plants <- function(n = 12L) {
  set.seed(10)
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

# --- ensureTxyColumns --------------------------------------------------------

test_that("ensureTxyColumns adds table_no/pos_x/pos_y", {
  raw <- data.frame(
    `T:X:Y` = c("1:2:3", "4:5:6"),
    unit = c("u1", "u2"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- ensureTxyColumns(raw)
  expect_true(all(c("table_no", "pos_x", "pos_y") %in% names(out)))
  expect_equal(out$table_no, c(1L, 4L))
  expect_equal(out$pos_x, c(2L, 5L))
  expect_equal(out$pos_y, c(3L, 6L))
})

# --- collapsePlantsBySelectedTimes -------------------------------------------

test_that("collapsePlantsBySelectedTimes yields 1 row/plant; selected ≠ all", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  plants <- length(unique(fix$unit))
  all_t <- sort(unique(fix$dbscan_cluster))
  selected <- c(1L, 2L)
  expect_lt(length(selected), length(all_t))

  subset <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = selected,
    metrics = c("mean", "auc", "peak")
  )
  full <- collapsePlantsBySelectedTimes(
    data = fix,
    trait = "height",
    time_col = "dbscan_cluster",
    time_values = all_t,
    metrics = c("mean", "auc", "peak")
  )

  expect_equal(nrow(subset), plants)
  expect_equal(length(unique(subset$unit)), plants)
  expect_equal(nrow(full), plants)
  expect_false(isTRUE(all.equal(subset$mean, full$mean)))
})

# --- residualizeSpatial ------------------------------------------------------

test_that("residualizeSpatial returns list(data, meta); passthrough without X/Y", {
  raw <- data.frame(
    height = c(1, 2, 3, 4),
    table_no = c(1L, 1L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  out <- residualizeSpatial(raw, y_col = "height")
  expect_type(out, "list")
  expect_true(all(c("data", "meta") %in% names(out)))
  expect_false(isTRUE(out$meta$spatial_applied))
  expect_equal(out$meta$resid_col, "height_spat_resid")
  expect_equal(out$data$height_spat_resid, raw$height)
})

test_that("residualizeSpatial with X/Y sets resid col when applied", {
  skip_if_not_installed("mgcv")
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  out <- residualizeSpatial(fix, y_col = "height")
  expect_type(out, "list")
  expect_true(all(c("data", "meta") %in% names(out)))
  expect_true(isTRUE(out$meta$spatial_applied))
  expect_true("height_spat_resid" %in% names(out$data))
  expect_equal(out$meta$resid_col, "height_spat_resid")
})

# --- fitBlupModels -----------------------------------------------------------

test_that("fitBlupModels return list has required names and geno columns", {
  skip_if_not_installed("lme4")
  collapsed <- .sf_iface_collapsed_mini()
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE
  )
  expect_true(all(c(
    "genotype_table", "treatment_table", "varcomp", "formula", "meta", "h2"
  ) %in% names(out)))
  expect_true(isTRUE(out$success))
  expect_true(all(c("level", "blup", "rank") %in% names(out$genotype_table)))
  if ("se" %in% names(out$genotype_table)) {
    expect_true(is.numeric(out$genotype_table$se))
  }
  expect_gt(nrow(out$genotype_table), 0L)
})

test_that("fitBlupModels multi-time guard fails with success=FALSE", {
  fix <- readRDS(sf_fixture_path("plant_curve_mini.rds"))
  out <- fitBlupModels(
    data = fix,
    trait = "height",
    apply_spatial = FALSE
  )
  expect_false(isTRUE(out$success))
  expect_true(is.character(out$reason) || is.character(out$meta$reason))
  expect_true(isTRUE(out$meta$n_time_per_plant > 1L))
  expect_equal(nrow(out$genotype_table), 0L)
})

test_that("fitBlupModels 1-level treatment yields empty treatment panel", {
  skip_if_not_installed("lme4")
  collapsed <- .sf_iface_collapsed_mini()
  expect_equal(length(unique(collapsed$treatment)), 1L)
  out <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE
  )
  expect_equal(nrow(out$treatment_table), 0L)
  expect_true(all(c("level", "estimate", "se", "rank") %in% names(out$treatment_table)))
  expect_match(out$meta$treatment_message, "treatment|level", ignore.case = TRUE)
})

test_that("fitBlupModels spatial passthrough when no X/Y keeps trait response", {
  skip_if_not_installed("lme4")
  raw <- data.frame(
    mean = c(10, 12, 14, 16, 11, 13, 15, 17),
    cultivar = rep(c("cvA", "cvB", "cvC", "cvD"), each = 2L),
    treatment = "ctrl",
    `repeat` = rep(1:2, 4L),
    table_no = rep(1:2, each = 4L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- fitBlupModels(
    data = raw,
    trait = "mean",
    apply_spatial = TRUE
  )
  expect_true(isTRUE(out$success))
  expect_false(isTRUE(out$meta$spatial_applied))
  expect_equal(out$meta$y_col, "mean")
  expect_match(out$formula, "^mean ~")
})

# --- reduceTraitsPCA ---------------------------------------------------------

test_that("reduceTraitsPCA: nrow scores == n plants; cultivar excluded; PCs bind", {
  plant_df <- .sf_iface_wide_plants()
  out <- reduceTraitsPCA(
    plant_df,
    trait_cols = c("height", "area", "volume", "ndvi", "cultivar"),
    n_comp = 2L
  )
  expect_true(isTRUE(out$meta$success))
  expect_equal(nrow(out$plant_df_with_pcs), nrow(plant_df))
  expect_false("cultivar" %in% out$meta$trait_cols)
  pc_cols <- grep("^PC[0-9]+$", names(out$plant_df_with_pcs), value = TRUE)
  expect_true(length(pc_cols) >= 1L)
  expect_true(all(c("unit", "cultivar", "PC1") %in% names(out$plant_df_with_pcs)))
  expect_equal(out$plant_df_with_pcs$unit, plant_df$unit)
})

# --- Committed interface contract fixture ------------------------------------

test_that("interface_contract.expected.rds freezes live formals and geometry_cols", {
  contract_path <- sf_fixture_path("interface_contract.expected.rds")
  expect_true(file.exists(contract_path))
  expected <- readRDS(contract_path)

  helper_fns <- list(
    ensureTxyColumns = ensureTxyColumns,
    collapsePlantsBySelectedTimes = collapsePlantsBySelectedTimes,
    residualizeSpatial = residualizeSpatial,
    fitBlupModels = fitBlupModels,
    reduceTraitsPCA = reduceTraitsPCA
  )
  expect_equal(sort(names(expected$helpers)), sort(names(helper_fns)))

  for (nm in names(helper_fns)) {
    live <- as.character(names(formals(helper_fns[[nm]])))
    frozen <- expected$helpers[[nm]]$formals
    expect_equal(live, frozen, info = paste0("formals drift: ", nm))
  }

  if (!is.null(expected$ensureTxyColumns$geometry_cols)) {
    expect_equal(
      expected$ensureTxyColumns$geometry_cols,
      c("table_no", "pos_x", "pos_y")
    )
  }
})

test_that("interface_contract.expected.rds matches live public shapes", {
  skip_if_not_installed("lme4")
  contract_path <- sf_fixture_path("interface_contract.expected.rds")
  expect_true(file.exists(contract_path))
  expected <- readRDS(contract_path)

  expect_true(all(c(
    "ensureTxyColumns",
    "collapsePlantsBySelectedTimes",
    "residualizeSpatial",
    "fitBlupModels",
    "reduceTraitsPCA"
  ) %in% names(expected$helpers)))

  collapsed <- .sf_iface_collapsed_mini()
  blup <- fitBlupModels(
    data = collapsed,
    trait = "mean",
    apply_spatial = FALSE
  )
  expect_equal(sort(names(blup)), sort(expected$fitBlupModels$return_names))
  expect_true(all(expected$fitBlupModels$genotype_cols %in% names(blup$genotype_table)))
  expect_equal(names(blup$treatment_table), expected$fitBlupModels$treatment_cols)

  spat <- residualizeSpatial(
    data.frame(height = 1:4, table_no = c(1L, 1L, 2L, 2L)),
    y_col = "height"
  )
  expect_equal(sort(names(spat)), sort(expected$residualizeSpatial$return_names))
  expect_true(all(expected$residualizeSpatial$meta_keys %in% names(spat$meta)))

  dr <- reduceTraitsPCA(.sf_iface_wide_plants(), n_comp = 2L)
  expect_equal(sort(names(dr)), sort(expected$reduceTraitsPCA$return_names))
})

# --- Optional soy-325 probe artifact (gitignored; skip if absent) ------------

test_that("optional soy325 probe RDS: success and ~325 genotypes", {
  probe_path <- file.path(sf_repo_root(), "tmp", "probe_blup_soy325.rds")
  skip_if_not(file.exists(probe_path), "tmp/probe_blup_soy325.rds not present")

  probe <- readRDS(probe_path)
  expect_true(isTRUE(probe$success))
  # top_k truncates genotype_table; full G count lives in meta
  n_g <- probe$meta$n_genotypes
  if (is.null(n_g)) {
    n_g <- nrow(probe$genotype_table)
  }
  expect_gte(as.integer(n_g), 300L)
  expect_lte(as.integer(n_g), 350L)
  expect_true(all(c("level", "blup", "rank") %in% names(probe$genotype_table)))
  expect_gt(nrow(probe$genotype_table), 0L)
})
