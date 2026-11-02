#!/usr/bin/env Rscript
# Biological narrative for soy-325 multi-trait PC BLUPs (+ single-time raw leaders).
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
})

source("shiny/common/pipeline_coach.R")
source("shiny/common/biology_readout.R")
source("shiny/common/trait_followup.R")
source("scripts/soy325_attach_cultivar_names.R")

out_dir <- "outputs/measurements/pc_blup_plots"
mt <- readRDS(file.path(out_dir, "multitrait_pc_blup.rds"))
project_id <- mt$project %||% "project_2025-07-14_b51-91_soy"
d <- as.data.frame(
  readRDS(file.path("data", project_id, paste0(project_id, "_merged_table.rds"))),
  stringsAsFactors = FALSE
)
d <- .sf_soy325_attach_cultivar_names(d, project_id)

tlev <- .sf_followup_time_levels(d)
time_slice <- if (length(tlev)) tlev[[length(tlev)]] else NULL

readout <- buildBiologyReadout(
  loadings = mt$loadings,
  variance = mt$variance,
  blup_by_pc = mt$blup_by_pc,
  check = mt$check,
  obs_data = d,
  top_k = 10L,
  overlap_k = 25L,
  time_value = time_slice
)

md_path <- file.path(out_dir, "BIOLOGY_READOUT_soy325.md")
writeLines(readout$markdown, md_path)
saveRDS(readout, file.path(out_dir, "biology_readout_soy325.rds"))

# Append short pointer on README
readme <- file.path(out_dir, "README_soy325.md")
if (file.exists(readme)) {
  writeLines(
    c(
      readLines(readme, warn = FALSE),
      "",
      "## Biology",
      "See `BIOLOGY_READOUT_soy325.md` — PC meaning, best/worst, cross-PC top-25, size↔NDVI groups."
    ),
    readme
  )
}

message("Wrote ", md_path)
cat("\n---- BIOLOGY BULLETS ----\n")
cat(paste0("- ", readout$bullets), sep = "\n")
cat("\n")
