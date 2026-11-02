#!/usr/bin/env Rscript
# Re-render honest multi-trait plots with readable labels (no overlap).
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(ggplot2)
})
source("shiny/common/plotting.R")
source("shiny/common/txy_parse.R")
source("shiny/common/plant_collapse.R")
source("shiny/common/trait_reduction.R")

out_dir <- "outputs/measurements/pc_blup_plots"
mt <- readRDS(file.path(out_dir, "multitrait_pc_blup.rds"))

save_p <- function(p, name, w, h) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 140, bg = "white", limitsize = FALSE)
  message(path)
}

save_p(
  createDrScreePlot(
    mt$variance,
    title = sprintf(
      "Scree — %d plant-mean phenotypes (PC1=%.0f%%)",
      mt$n_traits, 100 * mt$variance$prop_var[1]
    ),
    max_show = 12L
  ),
  "05_scree_multitrait.png", 9, 4.8
)

save_p(
  createDrLoadingsPlot(
    mt$loadings,
    max_pcs = 3L,
    top_n = 8L,
    title = "Loadings — multi-trait PCA (top-8 |loading| per PC)"
  ),
  "06_loadings_multitrait.png", 8.5, 14
)

d <- readRDS("data/project_2025-07-14_b51-91_soy/project_2025-07-14_b51-91_soy_merged_table.rds")
if (!"treatment" %in% names(d)) d$treatment <- "ctrl"
times <- sort(unique(d$dbscan_cluster))
sel <- times[c(1L, ceiling(length(times) / 2), length(times))]
sub <- d[d$dbscan_cluster %in% sel, , drop = FALSE]
plant <- sub |>
  dplyr::group_by(dplyr::across(dplyr::any_of(c("unit", "cultivar", "treatment")))) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(mt$trait_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )
dr <- reduceTraitsPCA(as.data.frame(plant), trait_cols = mt$trait_cols, scale = TRUE)
save_p(
  createDrScoresPlot(
    dr$plant_df_with_pcs,
    color_by = "cultivar",
    title = sprintf(
      "PC1 vs PC2 — %d traits (PC1=%.0f%%)",
      mt$n_traits, 100 * mt$variance$prop_var[1]
    )
  ),
  "07_scores_multitrait.png", 9, 6.5
)

n_geno <- nrow(mt$blup$genotype_table)
save_p(
  createBlupRankPlot(
    mt$blup$genotype_table,
    checks = mt$check,
    title = paste0("BLUP on multi-trait PC1 (check=", mt$check, ")")
  ),
  "08_blup_ranks_multitrait_pc1.png",
  8,
  max(5.5, 0.45 * n_geno + 2.2)
)

n_dun <- if (is.data.frame(mt$blup$meta$dunnett)) nrow(mt$blup$meta$dunnett) else 12L
save_p(
  createBlupDunnettPlot(
    mt$blup$meta$dunnett,
    title = "Dunnett vs check — multi-trait PC1"
  ),
  "09_dunnett_multitrait_pc1.png",
  8.5,
  max(5.5, 0.45 * n_dun + 2.4)
)

message("Readable multi-trait plots rewritten.")
