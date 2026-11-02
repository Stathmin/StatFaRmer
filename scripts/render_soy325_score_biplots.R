#!/usr/bin/env Rscript
# Re-render soy-325 PC score plots with biomass/NDVI GS-metric loading vectors.
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(ggplot2)
})
source("shiny/common/plotting.R")
source("shiny/common/txy_parse.R")
source("shiny/common/plant_collapse.R")
source("shiny/common/trait_reduction.R")
source("scripts/soy325_attach_cultivar_names.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b
out_dir <- "outputs/measurements/pc_blup_plots"
mt <- readRDS(file.path(out_dir, "multitrait_pc_blup.rds"))
project_id <- mt$project %||% "project_2025-07-14_b51-91_soy"

d <- as.data.frame(
  readRDS(file.path("data", project_id, paste0(project_id, "_merged_table.rds"))),
  stringsAsFactors = FALSE
)
if (!"treatment" %in% names(d)) d$treatment <- "ctrl"
d <- d[!is.na(d$cultivar), , drop = FALSE]
d <- .sf_soy325_attach_cultivar_names(d, project_id)
n_g <- length(unique(d$cultivar))

sel <- if (!is.null(mt$time_clusters) && length(mt$time_clusters)) {
  suppressWarnings(as.integer(mt$time_clusters))
} else {
  times <- sort(unique(d$dbscan_cluster))
  times[c(1L, ceiling(length(times) / 2), length(times))]
}
gs_metrics <- mt$gs_metrics %||% c("auc", "peak", "growth_rate", "last")
trait_pool <- mt$base_traits %||% {
  # fallback: traits embedded in metric column names
  unique(sub("__(auc|peak|growth_rate|last|mean)$", "", mt$trait_cols))
}
# Force biomass + NDVI into the pool (prior runs may have capped them out)
must <- intersect(c("digital_biomass_mm3", "ndvi_average"), names(d))
trait_pool <- unique(c(must, trait_pool))
if (length(trait_pool) > 12L) {
  rest <- setdiff(trait_pool, must)
  trait_pool <- c(must, rest)[seq_len(12L)]
}

plant <- collapseManyTraitsBySelectedTimes(
  data = d,
  traits = trait_pool,
  time_col = "dbscan_cluster",
  time_values = sel,
  metrics = gs_metrics
)
metric_cols <- attr(plant, "metric_cols")
dr_m <- reduceTraitsPCA(
  as.data.frame(plant),
  trait_cols = metric_cols,
  n_comp = NULL,
  scale = TRUE,
  require_named_temporal_metrics = TRUE
)
stopifnot(isTRUE(dr_m$meta$success))

pv <- dr_m$variance$prop_var
names(pv) <- as.character(dr_m$variance$component)
time_note <- sprintf(
  "GS metrics · clusters %s · not raw across-time means",
  paste(sel, collapse = ",")
)
message(time_note)

vec_traits <- c("digital_biomass_mm3", "ndvi_average")

save_p <- function(p, name, w = 10, h = 7.5) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 140, bg = "white", limitsize = FALSE)
  message("Wrote ", path)
}

for (pair in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))) {
  save_p(
    createDrScoresPlot(
      dr_m$plant_df_with_pcs,
      color_by = "cultivar",
      x_pc = pair[[1]],
      y_pc = pair[[2]],
      loadings = dr_m$loadings,
      vector_traits = vec_traits,
      subtitle = time_note,
      title = sprintf(
        "%s vs %s — %d GS metrics, %d genotypes (%.0f%% / %.0f%%)",
        pair[[1]], pair[[2]], length(metric_cols), n_g,
        100 * pv[[pair[[1]]]], 100 * pv[[pair[[2]]]]
      )
    ),
    sprintf("07_scores_%s_%s.png", tolower(pair[[1]]), tolower(pair[[2]]))
  )
}
file.copy(
  file.path(out_dir, "07_scores_pc1_pc2.png"),
  file.path(out_dir, "07_scores_multitrait.png"),
  overwrite = TRUE
)

mt$loadings <- dr_m$loadings
mt$variance <- dr_m$variance
mt$trait_cols <- dr_m$meta$trait_cols
mt$temporal_basis <- "named_GS_metrics"
mt$gs_metrics <- gs_metrics
mt$base_traits <- trait_pool
mt$time_clusters <- as.character(sel)
saveRDS(mt, file.path(out_dir, "multitrait_pc_blup.rds"))
message("Score biplots done (GS metrics).")
