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
  # Log resolved paths for diagnostics (use INFO so it appears by default)
  logEvent('INFO', 'minimal.generate.raw.paths', list(project_root = project_root, script = generate_raw_path))
  
  if (!file.exists(generate_raw_path)) {
    logError('generate_raw.R not found')
    stop('generate_raw.R not found')
  }
  
  # Create a temporary R script file to avoid shell escaping issues
  temp_script <- tempfile(fileext = '.R')
  script_content <- sprintf("
setwd('%s')
source('%s')
generateRawData('%s')
", project_root, generate_raw_path, project)
  writeLines(script_content, temp_script)
  logEvent('INFO', 'minimal.generate.raw.script', list(script_path = temp_script, content = script_content))
  
  # Preserve current library paths for the child so required packages are available without renv activation
  child_libs <- .libPaths()
  status <- tryCatch({
    system2('Rscript', args = temp_script, stdout = TRUE, stderr = TRUE, wd = project_root, env = c(paste0('R_LIBS=', paste(child_libs, collapse=':'))))
  }, error = function(e) {
    logError(paste('minimal.fallback.main.failed', e$message)); return(NA)
  }, finally = {
    unlink(temp_script)
  })
  exit_status <- tryCatch(attr(status, 'status'), error = function(e) NULL)
  if (length(status) == 0 || any(is.na(status))) {
    logError('minimal.generate.raw.failed.na_status')
  }
  logEvent('INFO', 'minimal.generate.raw.child_status', list(exit_status = ifelse(is.null(exit_status), 'NULL', as.character(exit_status))))
  if (length(status) > 0 && !any(is.na(status))) {
    logEvent('INFO', 'minimal.generate.raw.output', list(output = paste(status, collapse='\n')))
  }
  
  # Re-check project outputs
  raw_merged_path <- here::here('data', project, paste0(project, '_raw_merged_table.rds'))
  raw_groups_path <- here::here('data', project, paste0(project, '_raw_vector_of_groups.rds'))
  raw_merged_exists <- file.exists(raw_merged_path)
  raw_groups_exists <- file.exists(raw_groups_path)
  ok <- raw_merged_exists && raw_groups_exists
  if (!isTRUE(ok)) {
    logEvent('ERROR', 'minimal.generate.raw.outputs_missing', list(
      raw_merged = raw_merged_path,
      raw_groups = raw_groups_path,
      exists_merged = raw_merged_exists,
      exists_groups = raw_groups_exists
    ))
    # Fallback: try generating inline in current R session to avoid child env issues
    logEvent('WARN', 'minimal.generate.raw.inline_fallback.start', list(project = project))
    try({
      # Source into an isolated environment to avoid rm(list=ls()) nuking app symbols
      isolated_env <- new.env(parent = emptyenv())
      sys.source(generate_raw_path, envir = isolated_env)
      isolated_env$generateRawData(project)
    }, silent = TRUE)
    raw_merged_exists <- file.exists(raw_merged_path)
    raw_groups_exists <- file.exists(raw_groups_path)
    ok <- raw_merged_exists && raw_groups_exists
    logEvent(if (ok) 'INFO' else 'ERROR', 'minimal.generate.raw.inline_fallback.done', list(
      ok = ok,
      exists_merged = raw_merged_exists,
      exists_groups = raw_groups_exists
    ))
  }
  if (!isTRUE(ok)) stop('Minimal preprocess: generate_raw.R failed to create raw cache files')

  # Raw cache files are already in the correct location from generate_raw.R
  # No copying needed since generate_raw.R creates the _raw_ files directly

  # Inspect rows/cols
  raw_df <- tryCatch(readRDS(raw_merged), error = function(e) NULL)
  nrow_df <- if (!is.null(raw_df)) nrow(raw_df) else NA_integer_
  ncol_df <- if (!is.null(raw_df)) ncol(raw_df) else NA_integer_

  logEvent('INFO', 'minimal.raw.saved', list(project = project, raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
  return(list(raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
}


