#!/usr/bin/env Rscript
# Render PC/BLUP plots for chat review + optional multi-trait PCA contrast.
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(ggplot2)
})

source("shiny/common/plotting.R")
source("shiny/common/txy_parse.R")
source("shiny/common/plant_collapse.R")
source("shiny/common/spatial_residualize.R")
source("shiny/common/trait_reduction.R")
source("shiny/common/blup_ranking.R")
source("shiny/common/pipeline_coach.R")

out_dir <- "outputs/measurements/pc_blup_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_p <- function(p, name, w = 8, h = 5) {
  path <- file.path(out_dir, name)
  ggplot2::ggsave(path, p, width = w, height = h, dpi = 120, bg = "white")
  message("Wrote ", path)
  path
}

# Multi-trait plant PCA only. Prefer: Rscript scripts/render_soy325_pc_blup_suite.R
d <- readRDS("data/project_2025-07-14_b51-91_soy/project_2025-07-14_b51-91_soy_merged_table.rds")
if (!"treatment" %in% names(d)) d$treatment <- "ctrl"
times <- sort(unique(d$dbscan_cluster))
sel <- times[c(1L, ceiling(length(times) / 2), length(times))]

# Collapse several numeric canopy traits (not design cols) to plant means
design <- c(
  "unit", "V.T.R", "cultivar", "treatment", "repeat", "T:X:Y", "table_no",
  "pos_x", "pos_y", "dbscan_cluster", "timestamp_group", "outlier"
)
num <- names(d)[vapply(d, is.numeric, logical(1))]
trait_pool <- setdiff(num, design)
# drop near-constant / junk
keep <- vapply(trait_pool, function(nm) {
  x <- d[[nm]]
  is.numeric(x) && sum(is.finite(x)) > 50 && stats::sd(x, na.rm = TRUE) > 0
}, logical(1))
trait_pool <- trait_pool[keep]
# cap for speed/clarity
if (length(trait_pool) > 40L) {
  # prefer biomass/area/height-ish names first
  pref <- grep("biomass|area|height|volume|width|depth|leaf|ndvi|gli|hue|saturation",
               trait_pool, ignore.case = TRUE, value = TRUE)
  trait_pool <- unique(c(pref, trait_pool))[seq_len(40L)]
}

message("Multi-trait pool n=", length(trait_pool))

# One row per plant: mean over selected times for each trait
id <- if ("unit" %in% names(d)) "unit" else "V.T.R"
sub <- d[d$dbscan_cluster %in% sel, , drop = FALSE]
plant <- sub |>
  dplyr::group_by(dplyr::across(dplyr::any_of(c(id, "cultivar", "treatment", "repeat", "table_no", "pos_x", "pos_y")))) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trait_pool), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

dr_m <- reduceTraitsPCA(as.data.frame(plant), trait_cols = trait_pool, n_comp = NULL, scale = TRUE)
stopifnot(isTRUE(dr_m$meta$success))

p_scree_m <- createDrScreePlot(
  dr_m$variance,
  title = sprintf("Scree — plant means of %d phenotypes (honest DR)", length(dr_m$meta$trait_cols))
)
p_load_m <- createDrLoadingsPlot(
  dr_m$loadings,
  max_pcs = 3L,
  title = "Loadings — multi-trait plant PCA (PC1–PC3)"
)
p_scores_m <- createDrScoresPlot(
  dr_m$plant_df_with_pcs,
  color_by = "cultivar",
  title = sprintf(
    "PC1 vs PC2 — %d traits (PC1=%.0f%% var)",
    length(dr_m$meta$trait_cols),
    100 * dr_m$variance$prop_var[1]
  )
)

save_p(p_scree_m, "05_scree_multitrait.png")
save_p(p_load_m, "06_loadings_multitrait.png", w = 9, h = 7)
save_p(p_scores_m, "07_scores_multitrait.png", w = 8, h = 6)

# BLUP on multi-trait PC1
mns <- tapply(dr_m$plant_df_with_pcs[[dr_m$meta$trait_cols[1]]], dr_m$plant_df_with_pcs$cultivar, mean, na.rm = TRUE)
# mid check by PC1 itself
pc1_by_g <- tapply(dr_m$plant_df_with_pcs$PC1, dr_m$plant_df_with_pcs$cultivar, mean, na.rm = TRUE)
check <- names(sort(pc1_by_g))[ceiling(length(pc1_by_g) / 2)]
bl <- fitBlupModels(
  dr_m$plant_df_with_pcs,
  trait = "PC1",
  genotype_var = "cultivar",
  treatment_var = "treatment",
  apply_spatial = FALSE,
  checks = check,
  top_k = 40L
)
p_rank_m <- createBlupRankPlot(bl$genotype_table, checks = check,
                               title = paste0("BLUP on multi-trait PC1 (check=", check, ")"))
p_dun_m <- createBlupDunnettPlot(bl$meta$dunnett,
                                 title = "Dunnett vs check — multi-trait PC1")
save_p(p_rank_m, "08_blup_ranks_multitrait_pc1.png", h = 6)
save_p(p_dun_m, "09_dunnett_multitrait_pc1.png", h = 5.5)

saveRDS(
  list(
    n_traits = length(dr_m$meta$trait_cols),
    trait_cols = dr_m$meta$trait_cols,
    variance = dr_m$variance,
    loadings = dr_m$loadings,
    check = check,
    blup = bl,
    explain_pc1 = .sf_explain_pc(dr_m$loadings, "PC1", 5L),
    explain_pc2 = .sf_explain_pc(dr_m$loadings, "PC2", 5L)
  ),
  file.path(out_dir, "multitrait_pc_blup.rds")
)

message("Done. Plots in ", out_dir)
