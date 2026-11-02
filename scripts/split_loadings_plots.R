if (file.exists("renv/activate.R")) source("renv/activate.R")
suppressPackageStartupMessages(library(ggplot2))
source("shiny/common/plotting.R")

out <- "outputs/measurements/pc_blup_plots"
mt <- readRDS(file.path(out, "multitrait_pc_blup.rds"))

for (pc in c("PC1", "PC2", "PC3")) {
  ld2 <- mt$loadings[, c("trait", pc), drop = FALSE]
  p <- createDrLoadingsPlot(
    ld2,
    max_pcs = 1L,
    top_n = 12L,
    title = paste0("Loadings ", pc, " — multi-trait (top-12, readable)")
  )
  fn <- file.path(out, sprintf("06_loadings_%s.png", tolower(pc)))
  ggsave(fn, p, width = 9, height = 7.2, dpi = 140, bg = "white")
  message(fn)
}

ggsave(
  file.path(out, "05_scree_multitrait.png"),
  createDrScreePlot(
    mt$variance,
    title = sprintf("Scree — %d phenotypes (PC1=%.0f%%)", mt$n_traits, 100 * mt$variance$prop_var[1]),
    max_show = 12L
  ),
  width = 9, height = 4.8, dpi = 140, bg = "white"
)
message("done")
