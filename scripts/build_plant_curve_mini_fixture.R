# Regenerates tests/testthat/fixtures/plant_curve_mini.rds
# Tiny synthetic plant × time panel for BLUP/DR unit tests (no soy-325 RDS).
#
# Usage (repo root):
#   Rscript scripts/build_plant_curve_mini_fixture.R

root <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(root, "tests", "testthat", "fixtures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "plant_curve_mini.rds")

# 12 plants: 4 cultivars × 2 reps, spread across 3 tables with X/Y variation
plants <- data.frame(
  cultivar = c(
    "cvA", "cvB", "cvC", "cvD",
    "cvA", "cvB", "cvC", "cvD",
    "cvA", "cvB", "cvC", "cvD"
  ),
  treatment = "ctrl",
  `repeat` = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L),
  table_no = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L),
  pos_x = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 2L, 3L, 4L, 5L),
  pos_y = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

plants$`T:X:Y` <- paste(plants$table_no, plants$pos_x, plants$pos_y, sep = ":")
plants$unit <- plants$`T:X:Y`
plants$V.T.R <- paste(plants$cultivar, plants$treatment, plants$`repeat`, sep = ".")

# Three time clusters; deterministic growth curves per plant
clusters <- data.frame(
  dbscan_cluster = 1:3,
  timestamp_group = as.POSIXct(
    c("2025-01-01 10:00:00", "2025-01-08 10:00:00", "2025-01-15 10:00:00"),
    tz = "UTC"
  ),
  day = c(0, 7, 14),
  stringsAsFactors = FALSE
)

cultivar_base <- c(cvA = 10, cvB = 12, cvC = 14, cvD = 16)
rep_shift <- c(`1` = 0, `2` = 0.5)

rows <- list()
idx <- 0L
for (i in seq_len(nrow(plants))) {
  p <- plants[i, , drop = FALSE]
  base_h <- unname(cultivar_base[[p$cultivar]]) + unname(rep_shift[[as.character(p$`repeat`)]])
  base_a <- base_h * 1.5 + p$pos_x * 0.1 + p$pos_y * 0.05
  for (j in seq_len(nrow(clusters))) {
    idx <- idx + 1L
    t <- clusters[j, , drop = FALSE]
    height <- base_h + 0.8 * t$day + 0.05 * t$day^2 / 14
    area <- base_a + 1.2 * t$day + 0.02 * t$day * p$table_no
    rows[[idx]] <- data.frame(
      unit = p$unit,
      V.T.R = p$V.T.R,
      cultivar = p$cultivar,
      treatment = p$treatment,
      `repeat` = p$`repeat`,
      `T:X:Y` = p$`T:X:Y`,
      table_no = p$table_no,
      pos_x = p$pos_x,
      pos_y = p$pos_y,
      dbscan_cluster = t$dbscan_cluster,
      timestamp_group = t$timestamp_group,
      height = height,
      area = area,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

plant_curve_mini <- do.call(rbind, rows)
rownames(plant_curve_mini) <- NULL

saveRDS(plant_curve_mini, out_path)
message(
  sprintf(
    "Wrote %s (%d rows, %d plants, %d clusters)",
    out_path,
    nrow(plant_curve_mini),
    length(unique(plant_curve_mini$unit)),
    length(unique(plant_curve_mini$dbscan_cluster))
  )
)
