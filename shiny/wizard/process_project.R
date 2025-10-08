# Lightweight project processing for Shiny on-demand

library(here)

#' Ensure project RDS exists; if missing, run wizard/generate_raw.R for that project
#' @param project Project name (folder under data/)
ensureProjectRDS <- function(project) {
  merged_path <- here('data', project, paste0(project, '_merged_table.rds'))
  groups_path <- here('data', project, paste0(project, '_vector_of_groups.rds'))
  if (file.exists(merged_path) && file.exists(groups_path)) {
    return(invisible(TRUE))
  }
  # On first cache creation, do not generate here — let caller handle pipeline
  return(invisible(FALSE))
}

