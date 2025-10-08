# =============================================================================
# MINIMAL PREPROCESS (RAW CACHE)
# =============================================================================

#' Build minimal project cache (raw) using existing pipeline defaults
#' - Uses wizard/generate_raw.R with wizard-controlled steps disabled (no logit, no IQR)
#' - Writes raw cache under data/<project>/
#' @param project Project name (folder under data/)
#' @return Named list with paths to raw cache files
minimalPreprocess <- function(project) {
  # Ensure logging and utils are available when called outside app context
  try(suppressWarnings(source(here::here('shiny', 'common', 'logger.R'))), silent = TRUE)
  try(suppressWarnings(source(here::here('shiny', 'common', 'utils.R'))), silent = TRUE)

  raw_merged <- here::here('data', project, paste0(project, '_raw_merged_table.rds'))
  raw_groups <- here::here('data', project, paste0(project, '_raw_vector_of_groups.rds'))

  # If already present, skip
  if (file.exists(raw_merged) && file.exists(raw_groups)) {
    logEvent('INFO', 'minimal.raw.exists', list(project = project))
    return(list(raw_merged = raw_merged, raw_groups = raw_groups))
  }

  # Generate via wizard/generate_raw.R (replaces old main.R)
  logEvent('INFO', 'minimal.generate.raw', list(project = project))
  project_root <- here::here()
  generate_raw_path <- file.path(project_root, 'shiny', 'wizard', 'generate_raw.R')
  
  if (!file.exists(generate_raw_path)) {
    logError('generate_raw.R not found')
    stop('generate_raw.R not found')
  }
  
  # Create a temporary R script file to avoid shell escaping issues
  temp_script <- tempfile(fileext = '.R')
  script_content <- sprintf("
setwd('%s')
if (file.exists('%s')) source('%s')
Sys.setenv(PROJECT_NAME='%s')
source('%s')
", project_root, file.path(project_root, 'renv', 'activate.R'), file.path(project_root, 'renv', 'activate.R'), project, generate_raw_path)
  writeLines(script_content, temp_script)
  
  status <- tryCatch({
    system2('Rscript', args = temp_script, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    logError(paste('minimal.fallback.main.failed', e$message)); return(NA)
  }, finally = {
    unlink(temp_script)
  })
  if (length(status) == 0 || any(is.na(status))) {
    logError('minimal.generate.raw.failed.na_status')
  } else {
    logEvent('DEBUG', 'minimal.generate.raw.output', list(output = paste(status, collapse='\n')))
  }
  
  # Re-check project outputs
  ok <- file.exists(here::here('data', project, paste0(project, '_raw_merged_table.rds'))) &&
        file.exists(here::here('data', project, paste0(project, '_raw_vector_of_groups.rds')))
  if (!isTRUE(ok)) {
    stop('Minimal preprocess: generate_raw.R failed to create raw cache files')
  }

  # Raw cache files are already in the correct location from generate_raw.R
  # No copying needed since generate_raw.R creates the _raw_ files directly

  # Inspect rows/cols
  raw_df <- tryCatch(readRDS(raw_merged), error = function(e) NULL)
  nrow_df <- if (!is.null(raw_df)) nrow(raw_df) else NA_integer_
  ncol_df <- if (!is.null(raw_df)) ncol(raw_df) else NA_integer_

  logEvent('INFO', 'minimal.raw.saved', list(project = project, raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
  return(list(raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
}


