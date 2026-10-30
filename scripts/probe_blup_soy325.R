# Offline soy-325 BLUP probe — collapse all clusters, optional spatial, fitBlupModels.
# Writes gitignored artifact: tmp/probe_blup_soy325.rds
#
# Usage (repo root):
#   Rscript scripts/probe_blup_soy325.R
#   # or: "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" scripts/probe_blup_soy325.R
#
# Inspect (do not commit the RDS):
#   probe <- readRDS("tmp/probe_blup_soy325.rds")
#   str(probe, max.level = 2)
#   head(probe$genotype_table)

root <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

common <- file.path(root, "shiny", "common")
source(file.path(common, "txy_parse.R"), local = TRUE)
source(file.path(common, "plant_collapse.R"), local = TRUE)
source(file.path(common, "spatial_residualize.R"), local = TRUE)
source(file.path(common, "blup_ranking.R"), local = TRUE)

project_id <- "project_2025-07-14_b51-91_soy"
merged_path <- file.path(
  root, "data", project_id, paste0(project_id, "_merged_table.rds")
)
out_dir <- file.path(root, "tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "probe_blup_soy325.rds")

trait_raw <- "plant_height_averaged_mm"
metric <- "mean"
apply_spatial <- TRUE
top_k <- 20L

timings <- list()
errors <- character(0)
t0_all <- proc.time()[["elapsed"]]

probe_fail <- function(stage, err) {
  msg <- conditionMessage(err)
  errors <<- c(errors, sprintf("[%s] %s", stage, msg))
  message("ERROR ", stage, ": ", msg)
}

# --- load -------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
merged <- tryCatch(
  {
    raw <- readRDS(merged_path)
    as.data.frame(raw, stringsAsFactors = FALSE)
  },
  error = function(e) {
    probe_fail("load", e)
    NULL
  }
)
timings$load_s <- proc.time()[["elapsed"]] - t0

if (is.null(merged)) {
  saveRDS(
    list(
      success = FALSE,
      project_id = project_id,
      merged_path = merged_path,
      trait_raw = trait_raw,
      metric = metric,
      n_rows_merged = NA_integer_,
      n_plants = NA_integer_,
      n_clusters = NA_integer_,
      genotype_table = data.frame(),
      formula = NA_character_,
      h2 = NA_real_,
      fit = NULL,
      timings = timings,
      errors = errors
    ),
    out_path
  )
  quit(save = "no", status = 1L)
}

# Geometry: product ensureTxyColumns parses unit table:x:y when pos_x/pos_y absent.
merged <- tryCatch(
  ensureTxyColumns(merged),
  error = function(e) {
    probe_fail("ensureTxyColumns", e)
    merged
  }
)
if (!"V.T.R" %in% names(merged) && "v_t_r" %in% names(merged)) {
  merged[["V.T.R"]] <- as.character(merged$v_t_r)
}
if (!"treatment" %in% names(merged)) {
  # Soy-325 is single-treatment; placeholder keeps design col contract.
  merged$treatment <- "ctrl"
}

clusters <- sort(unique(merged$dbscan_cluster[!is.na(merged$dbscan_cluster)]))

# --- collapse (all clusters) ------------------------------------------------
t0 <- proc.time()[["elapsed"]]
plants <- tryCatch(
  collapsePlantsBySelectedTimes(
    data = merged,
    trait = trait_raw,
    time_col = "dbscan_cluster",
    time_values = clusters,
    metrics = metric
  ),
  error = function(e) {
    probe_fail("collapse", e)
    NULL
  }
)
timings$collapse_s <- proc.time()[["elapsed"]] - t0

# --- BLUP (spatial optional inside fitBlupModels) ---------------------------
fit <- NULL
t0 <- proc.time()[["elapsed"]]
if (!is.null(plants)) {
  fit <- tryCatch(
    fitBlupModels(
      data = plants,
      trait = metric,
      genotype_var = "cultivar",
      treatment_var = "treatment",
      repeat_var = "repeat",
      table_var = "table_no",
      apply_spatial = apply_spatial,
      checks = NULL,
      top_k = top_k
    ),
    error = function(e) {
      probe_fail("fitBlupModels", e)
      NULL
    }
  )
}
timings$fit_s <- proc.time()[["elapsed"]] - t0
timings$total_s <- proc.time()[["elapsed"]] - t0_all

geno <- if (!is.null(fit) && is.data.frame(fit$genotype_table)) {
  fit$genotype_table
} else {
  data.frame()
}

artifact <- list(
  success = isTRUE(fit$success),
  reason = if (!is.null(fit)) fit$reason else paste(errors, collapse = "; "),
  project_id = project_id,
  merged_path = merged_path,
  trait_raw = trait_raw,
  metric = metric,
  apply_spatial = apply_spatial,
  n_rows_merged = nrow(merged),
  n_plants = if (!is.null(plants)) nrow(plants) else NA_integer_,
  n_clusters = length(clusters),
  clusters = clusters,
  genotype_table = geno,
  genotype_head = utils::head(geno, 10L),
  n_genotypes = nrow(geno),
  formula = if (!is.null(fit)) fit$formula else NA_character_,
  h2 = if (!is.null(fit)) fit$h2 else NA_real_,
  varcomp = if (!is.null(fit)) fit$varcomp else NULL,
  meta = if (!is.null(fit)) fit$meta else NULL,
  fit = fit,
  timings = timings,
  errors = errors
)

saveRDS(artifact, out_path)

message(
  "Wrote ", out_path,
  " success=", artifact$success,
  " n_plants=", artifact$n_plants,
  " n_genotypes=", artifact$n_genotypes,
  " total_s=", round(timings$total_s, 2)
)
if (length(errors) > 0L) {
  message("Errors recorded: ", paste(errors, collapse = " | "))
}
