#!/usr/bin/env Rscript
# Add PC2/PC3 score + BLUP + Dunnett plots for soy-325 (reuse trait pool from suite RDS).
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(ggplot2)
})

source("shiny/common/plotting.R")
source("shiny/common/txy_parse.R")
source("shiny/common/plant_collapse.R")
source("shiny/common/trait_reduction.R")
source("shiny/common/blup_ranking.R")
source("shiny/common/pipeline_coach.R")
source("scripts/soy325_attach_cultivar_names.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- "outputs/measurements/pc_blup_plots"
project_id <- "project_2025-07-14_b51-91_soy"
mt_path <- file.path(out_dir, "multitrait_pc_blup.rds")
stopifnot(file.exists(mt_path))
mt <- readRDS(mt_path)

d <- as.data.frame(
  readRDS(file.path("data", project_id, paste0(project_id, "_merged_table.rds"))),
  stringsAsFactors = FALSE
)
if (!"treatment" %in% names(d)) d$treatment <- "ctrl"
d <- d[!is.na(d$cultivar), , drop = FALSE]
d <- .sf_soy325_attach_cultivar_names(d, project_id)
n_g <- length(unique(d$cultivar))
trait_pool <- mt$trait_cols

times <- sort(unique(d$dbscan_cluster))
sel <- times[c(1L, ceiling(length(times) / 2), length(times))]
id <- if ("unit" %in% names(d)) "unit" else "V.T.R"
sub <- d[d$dbscan_cluster %in% sel, , drop = FALSE]
plant <- sub |>
  dplyr::group_by(dplyr::across(dplyr::any_of(c(
    id, "cultivar", "treatment", "repeat", "table_no", "pos_x", "pos_y"
  )))) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trait_pool), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )
dr_m <- reduceTraitsPCA(as.data.frame(plant), trait_cols = trait_pool, n_comp = NULL, scale = TRUE)
stopifnot(isTRUE(dr_m$meta$success))
check_m <- mt$check

save_p <- function(p, name, w = 9, h = 6) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 140, bg = "white", limitsize = FALSE)
  message("Wrote ", path)
}

pv <- dr_m$variance$prop_var
names(pv) <- as.character(dr_m$variance$component)

time_note <- sprintf(
  "Plant means over clusters %s (NOT cluster-4-only)",
  paste(sel, collapse = ",")
)
for (pair in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))) {
  save_p(
    createDrScoresPlot(
      dr_m$plant_df_with_pcs,
      color_by = "cultivar",
      x_pc = pair[[1]],
      y_pc = pair[[2]],
      loadings = dr_m$loadings,
      vector_traits = c("digital_biomass_mm3", "ndvi_average"),
      subtitle = time_note,
      title = sprintf(
        "%s vs %s — %d traits, %d genotypes (%.0f%% / %.0f%%)",
        pair[[1]], pair[[2]], length(trait_pool), n_g,
        100 * pv[[pair[[1]]]], 100 * pv[[pair[[2]]]]
      )
    ),
    sprintf("07_scores_%s_%s.png", tolower(pair[[1]]), tolower(pair[[2]])),
    9, 6.5
  )
}
file.copy(
  file.path(out_dir, "07_scores_pc1_pc2.png"),
  file.path(out_dir, "07_scores_multitrait.png"),
  overwrite = TRUE
)

blup_by_pc <- list()
for (pc in c("PC1", "PC2", "PC3")) {
  message("BLUP ", pc, " …")
  bl <- fitBlupModels(
    dr_m$plant_df_with_pcs,
    trait = pc,
    genotype_var = "cultivar",
    treatment_var = "treatment",
    apply_spatial = FALSE,
    checks = check_m,
    top_k = 0L
  )
  blup_by_pc[[pc]] <- bl
  if (!isTRUE(bl$success)) {
    message("FAIL ", pc, ": ", bl$reason %||% "?")
    next
  }
  save_p(
    createBlupRankPlot(
      bl$genotype_table,
      checks = check_m,
      title = sprintf(
        "BLUP on multi-trait %s (%.0f%% var; check=%s)",
        pc, 100 * pv[[pc]], check_m
      ),
      keep_ends = 25L
    ),
    sprintf("08_blup_ranks_multitrait_%s.png", tolower(pc)),
    9, max(8, 0.22 * 50 + 3)
  )
  save_p(
    createBlupDunnettPlot(
      bl$meta$dunnett,
      title = sprintf("Dunnett vs check — multi-trait %s", pc),
      max_show = 50L,
      prefer_sig = TRUE
    ),
    sprintf("09_dunnett_multitrait_%s.png", tolower(pc)),
    9, 12
  )
}

mt$blup <- blup_by_pc$PC1
mt$blup_by_pc <- blup_by_pc
mt$explain_pc3 <- .sf_explain_pc(dr_m$loadings, "PC3", 5L)
saveRDS(mt, mt_path)

readme <- file.path(out_dir, "README_soy325.md")
if (file.exists(readme)) {
  extra <- c(
    "",
    "## PC2 / PC3",
    paste0("PC3: ", mt$explain_pc3),
    "- Scores: `07_scores_pc1_pc2.png`, `07_scores_pc1_pc3.png`, `07_scores_pc2_pc3.png`",
    "- BLUP ranks: `08_blup_ranks_multitrait_pc{1,2,3}.png`",
    "- Dunnett: `09_dunnett_multitrait_pc{1,2,3}.png`"
  )
  writeLines(c(readLines(readme, warn = FALSE), extra), readme)
}

message("PC2/PC3 plots done.")
