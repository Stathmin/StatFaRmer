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

  logEvent('INFO', 'minimal.generate.raw', list(project = project))
  project_root <- here::here()
  generate_raw_path <- file.path(project_root, 'shiny', 'wizard', 'generate_raw.R')
  logEvent('INFO', 'minimal.generate.raw.paths', list(project_root = project_root, script = generate_raw_path))

  if (!file.exists(generate_raw_path)) {
    logError('generate_raw.R not found')
    stop('generate_raw.R not found')
  }

  raw_merged_path <- raw_merged
  raw_groups_path <- raw_groups

  # --- 1) In-process first: same R / renv / .libPaths() as Shiny (avoids PATH Rscript mismatch). ---
  inline_err <- ''
  logEvent('INFO', 'minimal.generate.raw.inline_first', list(project = project))
  fb <- try({
    isolated_env <- new.env(parent = globalenv())
    sys.source(generate_raw_path, envir = isolated_env)
    isolated_env$generateRawData(project)
    invisible(NULL)
  }, silent = TRUE)
  if (inherits(fb, 'try-error')) {
    inline_err <- conditionMessage(attr(fb, 'condition'))
    logEvent('WARN', 'minimal.generate.raw.inline_first.failed', list(message = inline_err))
  }

  ok <- file.exists(raw_merged_path) && file.exists(raw_groups_path)

  # --- 2) Subprocess: Rscript from R.home() + explicit renv/activate.R (Rscript skips .Rprofile). ---
  child_log <- ''
  if (!isTRUE(ok)) {
    logEvent('INFO', 'minimal.generate.raw.subprocess', list(project = project))
    temp_script <- tempfile(fileext = '.R')
    script_root <- normalizePath(project_root, winslash = '/', mustWork = TRUE)
    script_gen <- normalizePath(generate_raw_path, winslash = '/', mustWork = TRUE)
    has_renv <- file.exists(file.path(project_root, 'renv', 'activate.R'))
    script_content <- paste0(
      'setwd(', encodeString(script_root, quote = '"'), ')\n',
      if (has_renv) 'source(file.path(getwd(), "renv", "activate.R"))\n' else '',
      'source(', encodeString(script_gen, quote = '"'), ')\n',
      'generateRawData(', encodeString(project, quote = '"'), ')\n',
      collapse = ''
    )
    writeLines(script_content, temp_script)
    logEvent('INFO', 'minimal.generate.raw.script', list(
      script_path = temp_script,
      renv_activate = has_renv
    ))

    rscript_cmd <- file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript')
    if (!file.exists(rscript_cmd)) {
      rscript_cmd <- 'Rscript'
    }

    child_libs <- .libPaths()
    status <- tryCatch({
      if (has_renv) {
        system2(rscript_cmd, args = temp_script, stdout = TRUE, stderr = TRUE, wd = project_root)
      } else {
        system2(
          rscript_cmd,
          args = temp_script,
          stdout = TRUE,
          stderr = TRUE,
          wd = project_root,
          env = c(paste0('R_LIBS=', paste(child_libs, collapse = .Platform$path.sep)))
        )
      }
    }, error = function(e) {
      logError(paste('minimal.generate.raw.subprocess.invoke_failed', e$message))
      NA
    }, finally = {
      unlink(temp_script)
    })

    exit_status <- tryCatch(attr(status, 'status'), error = function(e) NULL)
    if (length(status) > 0 && !any(is.na(status))) {
      child_log <- paste(status, collapse = '\n')
      logEvent('INFO', 'minimal.generate.raw.subprocess.output', list(output = child_log))
    } else {
      logError('minimal.generate.raw.subprocess.empty_or_na_status')
    }
    logEvent('INFO', 'minimal.generate.raw.child_status', list(
      exit_status = ifelse(is.null(exit_status), 'NULL', as.character(exit_status)),
      rscript_cmd = rscript_cmd
    ))

    ok <- file.exists(raw_merged_path) && file.exists(raw_groups_path)
    if (!isTRUE(ok) && nzchar(child_log)) {
      logEvent('ERROR', 'minimal.generate.raw.subprocess_failed_output', list(output = child_log))
    }
  }

  if (!isTRUE(ok)) {
    logEvent('ERROR', 'minimal.generate.raw.outputs_missing', list(
      raw_merged = raw_merged_path,
      raw_groups = raw_groups_path,
      exists_merged = file.exists(raw_merged_path),
      exists_groups = file.exists(raw_groups_path)
    ))
    detail <- trimws(paste(
      if (nzchar(inline_err)) paste0('In-session: ', inline_err),
      if (nzchar(child_log)) paste0('Rscript: ', substr(child_log, 1, 2000L)),
      sep = '\n'
    ))
    if (!nzchar(detail)) {
      detail <- 'No RDS written. Check logs/app.log for minimal.generate.raw.* events.'
    }
    stop(paste0('Minimal preprocess: generate_raw.R failed. ', detail))
  }

  raw_df <- tryCatch(readRDS(raw_merged), error = function(e) NULL)
  nrow_df <- if (!is.null(raw_df)) nrow(raw_df) else NA_integer_
  ncol_df <- if (!is.null(raw_df)) ncol(raw_df) else NA_integer_

  logEvent('INFO', 'minimal.raw.saved', list(project = project, raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
  return(list(raw_merged = raw_merged, raw_groups = raw_groups, rows = nrow_df, cols = ncol_df))
}
