#!/usr/bin/env Rscript
# Growth curves (full time) + single-time violins + Tukey letters.
# NEVER averages raw phenotypes across dates for the violin/Tukey panel.
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(ggplot2)
})

source("shiny/common/plotting.R")
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

out <- "outputs/measurements/pc_blup_plots"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Prefer full suite: Rscript scripts/render_soy325_pc_blup_suite.R
d <- readRDS("data/project_2025-07-14_b51-91_soy/project_2025-07-14_b51-91_soy_merged_table.rds")
d <- d[!is.na(d$cultivar), , drop = FALSE]

traits <- c(
  "digital_biomass_mm3",
  grep("^ndvi", names(d), value = TRUE)[1]
)
traits <- traits[!is.na(traits) & nzchar(traits) & traits %in% names(d)]

# Explicit last cluster — same default as UI; not a mean across clusters
tlev <- .sf_followup_time_levels(d)
time_slice <- if (length(tlev)) tlev[[length(tlev)]] else NULL
message("Violin/Tukey time slice: ", time_slice %||% "(none)")

save_p <- function(p, name, w = 10, h = 6) {
  path <- file.path(out, name)
  ggsave(path, p, width = w, height = h, dpi = 140, bg = "white", limitsize = FALSE)
  message(path)
}

for (tr in traits) {
  tag <- gsub("[^A-Za-z0-9]+", "_", tr)
  b <- buildTraitFollowupPlots(
    d, tr,
    group = "cultivar",
    time_value = time_slice,
    run_tukey = TRUE
  )
  message(tr, ": ", b$slice_note %||% "")
  save_p(b$curve, sprintf("10_curve_%s.png", tag), 11, 6.2)
  nlev <- length(unique(d$cultivar))
  save_p(b$violin, sprintf("11_violin_tukey_%s.png", tag), 11, max(6, 0.35 * nlev + 3))
  if (isTRUE(b$tukey$success)) {
    write.csv(
      b$tukey$letters,
      file.path(out, sprintf("11_letters_%s.csv", tag)),
      row.names = FALSE
    )
    message(
      "Tukey OK for ", tr, " @ ", b$time_value,
      " — ", nrow(b$tukey$letters), " letter rows; n plants=", nrow(b$tukey$plant_df)
    )
  } else {
    message("Tukey fail for ", tr, ": ", b$tukey$reason)
  }
}

message("Follow-up plots done (single-time violins only).")
