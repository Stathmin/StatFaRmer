#!/usr/bin/env Rscript
# Full PC / BLUP / follow-up suite for the hard case:
#   project_2025-07-14_b51-91_soy (~325 cultivars, no multi-treatment).
# Plots use high-cardinality presentation (ends / grey clouds / no Tukey letters).
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
source("shiny/common/letters_join.R")
source("shiny/common/letters_utils.R")
source("shiny/common/letters_normalize.R")
suppressWarnings({
  source("shiny/common/tukey_filter.R")
  source("shiny/common/letters_emmeans.R")
  source("shiny/common/letters_pvals.R")
  source("shiny/common/tukey.R")
})
source("shiny/common/trait_followup.R")
source("scripts/soy325_attach_cultivar_names.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

project_id <- "project_2025-07-14_b51-91_soy"
rds <- file.path("data", project_id, paste0(project_id, "_merged_table.rds"))
stopifnot(file.exists(rds))

out_dir <- "outputs/measurements/pc_blup_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/measurements", showWarnings = FALSE, recursive = TRUE)

save_p <- function(p, name, w = 9, h = 6) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 140, bg = "white", limitsize = FALSE)
  message("Wrote ", path)
  path
}

d <- readRDS(rds)
d <- as.data.frame(d, stringsAsFactors = FALSE)
if (!"treatment" %in% names(d)) d$treatment <- "ctrl"
d <- d[!is.na(d$cultivar), , drop = FALSE]
d <- .sf_soy325_attach_cultivar_names(d, project_id)
n_g <- length(unique(d$cultivar))
message(sprintf("Project %s: %d rows, %d cultivars, clusters=%s",
                project_id, nrow(d), n_g,
                paste(sort(unique(d$dbscan_cluster)), collapse = ",")))

# Time clusters for plant-level trait collapse (GS-style selected times)
times <- sort(unique(d$dbscan_cluster))
sel <- if (length(times) >= 3L) {
  times[c(1L, ceiling(length(times) / 2), length(times))]
} else {
  times
}

# ---------------------------------------------------------------------------
# Multi-trait plant PCA + BLUP on NAMED growth metrics (not across-time means)
# Literature: AUC / peak / rate / endpoint — never silent mean of raw scans.
# ---------------------------------------------------------------------------
design <- c(
  "unit", "V.T.R", "cultivar", "treatment", "repeat", "T:X:Y", "table_no",
  "pos_x", "pos_y", "dbscan_cluster", "timestamp_group", "outlier", "time_numeric"
)
num <- names(d)[vapply(d, is.numeric, logical(1))]
trait_pool <- setdiff(num, design)
keep <- vapply(trait_pool, function(nm) {
  x <- d[[nm]]
  is.numeric(x) && sum(is.finite(x)) > 50 && stats::sd(x, na.rm = TRUE) > 0
}, logical(1))
trait_pool <- trait_pool[keep]
# Cap base traits: each expands to 4 metrics (auc/peak/growth_rate/last)
pref <- grep(
  "biomass|area|height|volume|width|depth|leaf|ndvi|gli|hue|saturation|hull",
  trait_pool, ignore.case = TRUE, value = TRUE
)
# Always keep biomass + NDVI for biplot vectors / biology (must not be capped out)
must <- intersect(c("digital_biomass_mm3", "ndvi_average"), names(d))
trait_pool <- unique(c(must, pref, trait_pool))
if (length(trait_pool) > 12L) {
  rest <- setdiff(trait_pool, must)
  trait_pool <- c(must, rest)[seq_len(12L)]
}
gs_metrics <- c("auc", "peak", "growth_rate", "last")
message(
  "Multi-trait raw pool n=", length(trait_pool),
  " → metrics ", paste(gs_metrics, collapse = "/"),
  " over clusters ", paste(sel, collapse = ",")
)

plant <- collapseManyTraitsBySelectedTimes(
  data = d,
  traits = trait_pool,
  time_col = "dbscan_cluster",
  time_values = sel,
  metrics = gs_metrics
)
metric_cols <- attr(plant, "metric_cols")
.sf_dr_assert_named_temporal_metrics(metric_cols, context = "soy325 suite")
dr_m <- reduceTraitsPCA(
  as.data.frame(plant),
  trait_cols = metric_cols,
  n_comp = NULL,
  scale = TRUE,
  require_named_temporal_metrics = TRUE
)
stopifnot(isTRUE(dr_m$meta$success))

# Mid-ranked check on PC1 (shared across PC BLUPs for comparable Dunnett)
pc1_m <- tapply(dr_m$plant_df_with_pcs$PC1, dr_m$plant_df_with_pcs$cultivar, mean, na.rm = TRUE)
check_m <- names(sort(pc1_m))[ceiling(length(pc1_m) / 2)]

pcs_fit <- intersect(c("PC1", "PC2", "PC3"), names(dr_m$plant_df_with_pcs))
blup_by_pc <- list()
for (pc in pcs_fit) {
  message("BLUP multi-trait ", pc, " vs check=", check_m, " …")
  blup_by_pc[[pc]] <- fitBlupModels(
    dr_m$plant_df_with_pcs,
    trait = pc,
    genotype_var = "cultivar",
    treatment_var = "treatment",
    apply_spatial = FALSE,
    checks = check_m,
    top_k = 0L
  )
}
bl_m <- blup_by_pc[["PC1"]]

save_p(
  createDrScreePlot(
    dr_m$variance,
    title = sprintf(
      "Scree — %d GS metrics (%d traits × %s; PC1=%.0f%%)",
      length(dr_m$meta$trait_cols), length(trait_pool),
      paste(gs_metrics, collapse = "+"),
      100 * dr_m$variance$prop_var[1]
    ),
    max_show = 12L
  ),
  "05_scree_multitrait.png", 9, 4.8
)
save_p(
  createDrLoadingsPlot(
    dr_m$loadings, max_pcs = 3L, top_n = 8L,
    title = "Loadings — multi-trait GS-metric PCA (top-8 |loading| per PC)"
  ),
  "06_loadings_multitrait.png", 11, 14
)
for (pc in c("PC1", "PC2", "PC3")) {
  if (!pc %in% names(dr_m$loadings)) next
  ld2 <- dr_m$loadings[, c("trait", pc), drop = FALSE]
  save_p(
    createDrLoadingsPlot(ld2, max_pcs = 1L, top_n = 12L,
                         title = paste0("Loadings ", pc, " — multi-trait (top-12)")),
    sprintf("06_loadings_%s.png", tolower(pc)), 11, 7.5
  )
}

# Score scatters: PC1–PC2, PC1–PC3, PC2–PC3
score_pairs <- list(
  c("PC1", "PC2"),
  c("PC1", "PC3"),
  c("PC2", "PC3")
)
for (pair in score_pairs) {
  if (!all(pair %in% names(dr_m$plant_df_with_pcs))) next
  pv <- dr_m$variance$prop_var
  names(pv) <- as.character(dr_m$variance$component)
  # One arrow each: resolver maps base names → __auc (preferred)
  vec_traits <- c("digital_biomass_mm3", "ndvi_average")
  save_p(
    createDrScoresPlot(
      dr_m$plant_df_with_pcs,
      color_by = "cultivar",
      x_pc = pair[[1]],
      y_pc = pair[[2]],
      loadings = dr_m$loadings,
      vector_traits = vec_traits,
      subtitle = sprintf(
        "GS metrics (auc+peak+growth_rate+last) · clusters %s · not raw across-time means",
        paste(sel, collapse = ",")
      ),
      title = sprintf(
        "%s vs %s — %d GS metrics, %d genotypes (%.0f%% / %.0f%%)",
        pair[[1]], pair[[2]], length(dr_m$meta$trait_cols), n_g,
        100 * pv[[pair[[1]]]], 100 * pv[[pair[[2]]]]
      )
    ),
    sprintf("07_scores_%s_%s.png", tolower(pair[[1]]), tolower(pair[[2]])),
    9.5, 7.2
  )
}
# Keep legacy alias for PC1–PC2
if (file.exists(file.path(out_dir, "07_scores_pc1_pc2.png"))) {
  file.copy(
    file.path(out_dir, "07_scores_pc1_pc2.png"),
    file.path(out_dir, "07_scores_multitrait.png"),
    overwrite = TRUE
  )
}

for (pc in pcs_fit) {
  bl <- blup_by_pc[[pc]]
  if (!isTRUE(bl$success)) {
    message("Skip plots for ", pc, ": ", bl$reason %||% "BLUP failed")
    next
  }
  pc_l <- tolower(pc)
  pv <- dr_m$variance$prop_var[match(pc, dr_m$variance$component)]
  save_p(
    createBlupRankPlot(
      bl$genotype_table, checks = check_m,
      title = sprintf("BLUP on multi-trait %s (%.0f%% var; check=%s)", pc, 100 * pv, check_m),
      keep_ends = 25L
    ),
    sprintf("08_blup_ranks_multitrait_%s.png", pc_l),
    12, max(9, 0.24 * 50 + 3)
  )
  dun_ok <- tryCatch({
    save_p(
      createBlupDunnettPlot(
        bl$meta$dunnett,
        title = sprintf("Dunnett vs check — multi-trait %s", pc),
        max_show = 50L, prefer_sig = TRUE
      ),
      sprintf("09_dunnett_multitrait_%s.png", pc_l),
      12, 13
    )
    TRUE
  }, error = function(e) {
    message("Dunnett plot skipped for ", pc, ": ", conditionMessage(e))
    FALSE
  })
  invisible(dun_ok)
}

mt_bundle <- list(
  project = project_id,
  n_genotypes = n_g,
  n_traits = length(dr_m$meta$trait_cols),
  n_base_traits = length(trait_pool),
  base_traits = trait_pool,
  gs_metrics = gs_metrics,
  trait_cols = dr_m$meta$trait_cols,
  time_clusters = as.character(sel),
  temporal_basis = "named_GS_metrics",
  variance = dr_m$variance,
  loadings = dr_m$loadings,
  check = check_m,
  blup = bl_m,
  blup_by_pc = blup_by_pc,
  explain_pc1 = .sf_explain_pc(dr_m$loadings, "PC1", 5L),
  explain_pc2 = .sf_explain_pc(dr_m$loadings, "PC2", 5L),
  explain_pc3 = .sf_explain_pc(dr_m$loadings, "PC3", 5L)
)
saveRDS(mt_bundle, file.path(out_dir, "multitrait_pc_blup.rds"))

sig_from_bl <- function(bl) {
  dun <- bl$meta$dunnett
  if (!is.data.frame(dun) || nrow(dun) < 1L) return(dun)
  pcol <- intersect(c("p.value", "p.adj"), names(dun))[1L]
  if (is.na(pcol)) return(dun)
  dun[is.finite(dun[[pcol]]) & dun[[pcol]] < 0.05, , drop = FALSE]
}
sig_m <- sig_from_bl(bl_m)

saveRDS(
  list(
    project = project_id,
    n_traits = mt_bundle$n_traits,
    time_clusters = as.character(sel),
    n_plants = nrow(dr_m$plant_df_with_pcs),
    n_genotypes = n_g,
    variance = dr_m$variance,
    loadings = dr_m$loadings,
    pc1_meaning = mt_bundle$explain_pc1,
    pc2_meaning = mt_bundle$explain_pc2,
    pc3_meaning = mt_bundle$explain_pc3,
    check = check_m,
    blup_pc1_ok = isTRUE(bl_m$success),
    geno_pc1 = bl_m$genotype_table,
    dunnett_pc1 = bl_m$meta$dunnett,
    sig_pc1 = sig_m,
    geno_pc2 = if (isTRUE(blup_by_pc$PC2$success)) blup_by_pc$PC2$genotype_table else NULL,
    dunnett_pc2 = if (isTRUE(blup_by_pc$PC2$success)) blup_by_pc$PC2$meta$dunnett else NULL,
    sig_pc2 = if (isTRUE(blup_by_pc$PC2$success)) sig_from_bl(blup_by_pc$PC2) else NULL,
    geno_pc3 = if (isTRUE(blup_by_pc$PC3$success)) blup_by_pc$PC3$genotype_table else NULL,
    dunnett_pc3 = if (isTRUE(blup_by_pc$PC3$success)) blup_by_pc$PC3$meta$dunnett else NULL,
    sig_pc3 = if (isTRUE(blup_by_pc$PC3$success)) sig_from_bl(blup_by_pc$PC3) else NULL
  ),
  "outputs/measurements/pc_blup_demo_results.rds"
)

# ---------------------------------------------------------------------------
# Follow-up: curves + single-time ranked strips (no across-time mean; no Tukey)
# ---------------------------------------------------------------------------
tlev <- .sf_followup_time_levels(d)
time_slice <- if (length(tlev)) tlev[[length(tlev)]] else NULL
traits_fu <- c("digital_biomass_mm3", grep("^ndvi_average$", names(d), value = TRUE)[1])
traits_fu <- traits_fu[!is.na(traits_fu) & nzchar(traits_fu) & traits_fu %in% names(d)]

for (tr in traits_fu) {
  tag <- gsub("[^A-Za-z0-9]+", "_", tr)
  b <- buildTraitFollowupPlots(
    d, tr, group = "cultivar", time_value = time_slice, run_tukey = TRUE
  )
  message(tr, ": ", b$slice_note %||% "", " | Tukey: ", b$tukey$reason %||% "ok")
  save_p(b$curve, sprintf("10_curve_%s.png", tag), 11, 6.2)
  # ends plot is tall (40+40); wider for Cyrillic cultivar_name labels
  save_p(b$violin, sprintf("11_violin_tukey_%s.png", tag), 12, 16)
  if (isTRUE(b$tukey$success) && !is.null(b$tukey$letters)) {
    write.csv(b$tukey$letters, file.path(out_dir, sprintf("11_letters_%s.csv", tag)), row.names = FALSE)
  } else {
    writeLines(
      c(
        sprintf("project: %s", project_id),
        sprintf("trait: %s", tr),
        sprintf("time_slice: %s", time_slice %||% "NA"),
        sprintf("n_genotypes: %d", n_g),
        sprintf("tukey: %s", b$tukey$reason %||% "n/a"),
        "Presentation: ranked top/bottom genotype ends (createTraitRankStripPlot)."
      ),
      file.path(out_dir, sprintf("11_letters_%s_SKIPPED.txt", tag))
    )
  }
}

# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
writeLines(c(
  paste("# PC/BLUP plots —", project_id),
  "",
  sprintf("- **Genotypes:** %d (hard screen; no multi-level treatment).", n_g),
  sprintf("- **Time clusters used for GS collapse:** %s", paste(sel, collapse = ", ")),
  sprintf("- **Follow-up violin time slice:** %s (single-time; not across-time mean).", time_slice %||% "n/a"),
  "",
  "## Presentation rules for many genotypes",
  "- Growth curves: 10–90% ribbon + red median + sample trajectories (no 325-color legend).",
  "- Single-time distributions: top/bottom 40 genotype means (horizontal boxes), not 325 violins.",
  "- Tukey CLD: skipped when genotype levels > 20; use BLUP ranks + Dunnett vs check.",
  "- BLUP forests / Dunnett: top/bottom ends or sig + largest |Δ| (max ~50 labels).",
  "- PC scores: grey cloud (cultivar color legend suppressed).",
  "",
  sprintf(
    "Multi-trait PCA: **%d** named GS metrics from %d traits (%s) over clusters %s; PC1=%.1f%%; PC2=%.1f%%.",
    length(dr_m$meta$trait_cols), length(trait_pool),
    paste(gs_metrics, collapse = "+"),
    paste(sel, collapse = ","),
    100 * dr_m$variance$prop_var[1],
    100 * dr_m$variance$prop_var[2]
  ),
  paste0("PC1: ", .sf_explain_pc(dr_m$loadings, "PC1", 5L)),
  paste0("PC2: ", .sf_explain_pc(dr_m$loadings, "PC2", 5L)),
  sprintf("Dunnett check: %s; BLUP success: %s", check_m, isTRUE(bl_m$success)),
  "",
  "See `README_temporal_DR_basis.md` for literature basis (AUC/rate/endpoint, not raw means)."
), file.path(out_dir, "README_soy325.md"))

md <- "outputs/measurements/pc_blup_demo_results.md"
writeLines(c(
  "# Multi-trait PC BLUP (soy-325)",
  "",
  sprintf("- **Project:** %s", project_id),
  sprintf("- **Traits in PCA:** %d", length(dr_m$meta$trait_cols)),
  sprintf("- **Time clusters:** %s", paste(sel, collapse = ", ")),
  sprintf("- **Plants / genotypes:** %d / %d", nrow(dr_m$plant_df_with_pcs), n_g),
  sprintf("- **PC1 / PC2:** %.1f%% / %.1f%%", 100 * dr_m$variance$prop_var[1], 100 * dr_m$variance$prop_var[2]),
  sprintf("- **Dunnett check:** %s", check_m),
  "",
  "## PC meaning",
  sprintf("- **PC1:** %s", mt_bundle$explain_pc1),
  sprintf("- **PC2:** %s", mt_bundle$explain_pc2),
  "",
  sprintf("Significant Dunnett (p < 0.05) on PC1: %s",
          if (is.data.frame(sig_m) && nrow(sig_m) > 0) paste(nrow(sig_m), "contrasts") else "none / unavailable")
), md)

# Drop leftover pedagogical 3-metric plots if present
old_junk <- file.path(out_dir, c(
  "01_scree_3metrics.png", "02_loadings_3metrics.png",
  "03_blup_ranks_pc1_3metrics.png", "04_dunnett_pc1_3metrics.png",
  "README_why_99pct_was_wrong.md"
))
unlink(old_junk[file.exists(old_junk)])

# Biological narrative (axes, best/worst, cross-PC, size↔NDVI groups)
if (file.exists("shiny/common/biology_readout.R")) {
  source("shiny/common/biology_readout.R")
  source("shiny/common/trait_followup.R")
  bio <- buildBiologyReadout(
    loadings = dr_m$loadings,
    variance = dr_m$variance,
    blup_by_pc = blup_by_pc,
    check = check_m,
    obs_data = d,
    top_k = 10L,
    overlap_k = 25L,
    time_value = .sf_followup_time_levels(d) |> (\(x) if (length(x)) x[[length(x)]] else NULL)()
  )
  writeLines(bio$markdown, file.path(out_dir, "BIOLOGY_READOUT_soy325.md"))
  saveRDS(bio, file.path(out_dir, "biology_readout_soy325.rds"))
}

message("Soy-325 suite done → ", out_dir)
