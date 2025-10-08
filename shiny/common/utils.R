# StatFaRmer Utility Functions
# Iteration 3: Modular Shiny

library(here)
library(dbscan)
suppressed <- try({
  suppressPackageStartupMessages(library(jsonlite))
}, silent = TRUE)

# =============================================================================
# PROJECT MANAGEMENT
# =============================================================================

#' Get available projects with deploy restrictions
#' @return Vector of available project names
getAvailableProjects <- function() {
  # Load deploy config if it exists
  deploy_config_path <- here('config', 'deploy_config.R')
  if (file.exists(deploy_config_path)) {
    source(deploy_config_path)
    # Filter projects by allowed_projects if in deploy mode
    if (exists('allowed_projects')) {
      all_projects <- list.dirs(here('data'), recursive = FALSE, full.names = FALSE)
      all_projects <- all_projects[!grepl('^\\.', all_projects)]
      return(intersect(all_projects, allowed_projects))
    }
  }
  
  # Return all projects if no deploy restrictions
  all_projects <- list.dirs(here('data'), recursive = FALSE, full.names = FALSE)
  return(all_projects[!grepl('^\\.', all_projects)])
}

#' Select default project honoring environment override
#' @param available_projects character vector of project names
#' @return Single project name
getDefaultProject <- function(available_projects) {
  if (length(available_projects) == 0) return(NA_character_)
  selected <- Sys.getenv('SELECTED_PROJECT', '')
  if (nzchar(selected) && selected %in% available_projects) {
    return(selected)
  }
  if ("project_NO3" %in% available_projects) {
    return("project_NO3")
  }
  return(available_projects[1])
}

# =============================================================================
# DEBUG CACHE (24h TTL)
# =============================================================================

#' Get cache directory path under logs/cache
#' @return string path
getCacheDir <- function() {
  here('logs', 'cache')
}

#' Build cache file path for a given key
#' @param key unique cache key (safe for filenames)
#' @return absolute file path
cachePathFor <- function(key) {
  dir.create(getCacheDir(), showWarnings = FALSE, recursive = TRUE)
  here(getCacheDir(), paste0(key, '.rds'))
}

#' Save an object to cache when ENABLE_DEBUG_CACHE is TRUE
#' @param key cache key
#' @param object any R object
saveDebugCache <- function(key, object) {
  enable <- tryCatch(get('ENABLE_DEBUG_CACHE', envir = .GlobalEnv), error = function(e) FALSE)
  if (!isTRUE(enable)) return(invisible(FALSE))
  path <- cachePathFor(key)
  tryCatch({
    saveRDS(object, path)
    TRUE
  }, error = function(e) FALSE)
}

#' Load an object from cache if present and fresh (≤24h)
#' @param key cache key
#' @param ttl_hours numeric TTL hours (default 24)
#' @return cached object or NULL
loadDebugCache <- function(key, ttl_hours = 24) {
  enable <- tryCatch(get('ENABLE_DEBUG_CACHE', envir = .GlobalEnv), error = function(e) FALSE)
  if (!isTRUE(enable)) return(NULL)
  path <- cachePathFor(key)
  if (!file.exists(path)) return(NULL)
  info <- file.info(path)
  age_hours <- as.numeric(difftime(Sys.time(), info$mtime, units = 'hours'))
  if (is.na(age_hours) || age_hours > ttl_hours) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

#' Load project data with fallback (wizard helper)
#' 
#' Unified function for wizard data loading with raw/processed fallback.
#' Tries raw data first (for wizard), falls back to processed if needed.
#' 
#' @param project_name Name of the project to load
#' @param prefer 'raw' or 'processed' - which to try first
#' @param context Context string for logging (e.g., 'plotData')
#' @return Data frame or NULL if failed
loadProjectDataWithFallback <- function(project_name, prefer = 'raw', context = 'data') {
  if (!nzchar(project_name)) {
    logError(paste(context, ': empty project name'))
    return(NULL)
  }
  
  # Determine paths
  raw_path <- here('data', project_name, paste0(project_name, '_raw_merged_table.rds'))
  processed_path <- here('data', project_name, paste0(project_name, '_merged_table.rds'))
  
  # Try preferred path first
  first_path <- if (prefer == 'raw') raw_path else processed_path
  second_path <- if (prefer == 'raw') processed_path else raw_path
  
  # Try first choice
  if (file.exists(first_path)) {
    df <- tryCatch(readRDS(first_path), error = function(e) {
      logError(paste(context, '.read_failed:', e$message))
      NULL
    })
    if (!is.null(df)) {
      logEvent('DEBUG', paste(context, '.loaded'), list(
        project = project_name,
        source = prefer,
        nrow = nrow(df),
        ncol = ncol(df)
      ))
      return(df)
    }
  }
  
  # Fallback to second choice
  if (file.exists(second_path)) {
    df <- tryCatch(readRDS(second_path), error = function(e) {
      logError(paste(context, '.fallback_failed:', e$message))
      NULL
    })
    if (!is.null(df)) {
      logEvent('DEBUG', paste(context, '.loaded_fallback'), list(
        project = project_name,
        source = ifelse(prefer == 'raw', 'processed', 'raw'),
        nrow = nrow(df),
        ncol = ncol(df)
      ))
      return(df)
    }
  }
  
  # Both failed
  logError(paste(context, ': No raw or processed data available for', project_name))
  return(NULL)
}

#' Load and prepare data for a specific project
#' @param project_name Name of the project to load
#' @return List with merged_table and vector_of_groups
loadProjectData <- function(project_name) {
  
  # Check if project exists
  available_projects <- getAvailableProjects()
  if (!project_name %in% available_projects) {
    stop(paste("Project", project_name, "not available or not allowed"))
  }
  
  # Check if .rds files exist for this project (project folder only)
  merged_table_path <- here('data', project_name, paste0(project_name, '_merged_table.rds'))
  vector_of_groups_path <- here('data', project_name, paste0(project_name, '_vector_of_groups.rds'))
  
  if (!(file.exists(merged_table_path) && file.exists(vector_of_groups_path))) {
    stop(paste(
      "Data files not found for project", project_name, 
      "- expected:", basename(merged_table_path), "and", basename(vector_of_groups_path),
      "in", here('data', project_name),
      "\nPlease run Master Wizard."
    ))
  }
  
  # Load data from .rds files
  merged_table <- readRDS(merged_table_path)
  vector_of_groups <- readRDS(vector_of_groups_path)
  
  # Clean infinite values
  merged_table[sapply(merged_table, is.infinite)] <- NA

  # Convert metadata columns to factors (convention: all metadata as factors)
  metadata_cols <- c('timestamp', 'unit', 'v_t_r', 'treatment', 'cultivar', 'dbscan_cluster', 'timestamp_group')
  for (col in metadata_cols) {
    if (col %in% names(merged_table) && is.character(merged_table[[col]])) {
      merged_table[[col]] <- as.factor(merged_table[[col]])
    }
  }

  # Legacy compatibility: provide timestamp_group if only dbscan_cluster exists
  if (!("timestamp_group" %in% names(merged_table)) && ("dbscan_cluster" %in% names(merged_table))) {
    merged_table$timestamp_group <- as.factor(as.character(merged_table$dbscan_cluster))
  }
  
  return(list(
    merged_table = merged_table,
    vector_of_groups = vector_of_groups
  ))
}

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

#' Load and prepare data for Shiny application
#' @return List with merged_table and vector_of_groups
loadShinyData <- function() {
  
  # Load from .rds files (created by main.R)
  merged_table_path <- here('shiny', 'merged_table.rds')
  vector_of_groups_path <- here('shiny', 'vector_of_groups.rds')
  
  # Check if .rds files exist
  if (!file.exists(merged_table_path)) {
    stop("merged_table.rds not found. Please run src/main.R first to create data files.")
  }
  if (!file.exists(vector_of_groups_path)) {
    stop("vector_of_groups.rds not found. Please run src/main.R first to create data files.")
  }
  
  # Load data from .rds files
  merged_table <- readRDS(merged_table_path)
  vector_of_groups <- readRDS(vector_of_groups_path)
  
  # Clean infinite values
  merged_table[sapply(merged_table, is.infinite)] <- NA

  # Legacy compatibility: provide timestamp_group if only dbscan_cluster exists
  if (!("timestamp_group" %in% names(merged_table)) && ("dbscan_cluster" %in% names(merged_table))) {
    merged_table$timestamp_group <- as.character(merged_table$dbscan_cluster)
  }
  
  return(list(
    merged_table = merged_table,
    vector_of_groups = vector_of_groups
  ))
}

#' Get unique values from a table column
#' @param table Data frame
#' @param string Column name
#' @return Sorted unique values
getUnique <- function(table, string) {
  if (!string %in% names(table)) return(character(0))
  table %>%
    dplyr::select(dplyr::all_of(string)) %>%
    dplyr::distinct() %>%
    dplyr::pull(1) %>%
    sort()
}

#' Create ANOVA formula
#' @param response_var Response variable
#' @param factors ANOVA factors
#' @return Formula string
formulate <- function(response_var, factors) {
  if (length(factors) == 0) {
    return(paste(response_var, "~ 1"))
  }
  
  formula_str <- paste(response_var, "~", paste(factors, collapse = " + "))
  return(formula_str)
}

# =============================================================================
# STATISTICAL FUNCTIONS
# =============================================================================

#' Create ANOVA formula
#' @param RHS Response variable
#' @param LHF Left-hand factors
#' @return Formula string
formulate <- function(RHS, LHF) {
  if (length(LHF) == 0) {
    LHS = '1'
  } else if ((length(LHF) == 1)) {
    LHS = str_interp('1 + ${LHF}')
  } else {
    LHS = str_interp('1 + (${paste(LHF, collapse=" + ")})^2')
  }
  return(paste(RHS, LHS, sep = " ~ "))
}

#' Generate Tukey HSD labels
#' @param TUKEY Tukey test results
#' @param variable Variable name
#' @return Data frame with labels
generateLabelDf <- function(TUKEY, variable) {
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][, 4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Keep labels in the same order as in the boxplot
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment), ]
  colnames(Tukey.labels) = c('letter', 'group')
  return(Tukey.labels)
}

# =============================================================================
# DATA PREPARATION FOR UI
# =============================================================================

#' Prepare data selections for UI
#' @param merged_table Main data table
#' @return List with UI data selections
prepareUISelections <- function(merged_table) {
  # ANOVA factors - get character columns (excluding timestamp); drop 'unit'
  char_cols <- names(merged_table)[sapply(merged_table, function(x) is.character(x) || is.factor(x))]
  anova_factors <- sort(setdiff(char_cols, 'unit'))
  # Ensure dbscan_cluster is available as a factor option
  if ('dbscan_cluster' %in% names(merged_table)) {
    anova_factors <- sort(unique(c(anova_factors, 'dbscan_cluster')))
  }
  
  # Basic selections with schema fallbacks
  treatments <- getUnique(merged_table, 'treatment')
  cultivar_col <- if ('cultivar' %in% names(merged_table)) 'cultivar' else if ('genotype' %in% names(merged_table)) 'genotype' else NULL
  cultivars <- if (!is.null(cultivar_col)) getUnique(merged_table, cultivar_col) else character(0)
  
  # Use dbscan_cluster as timestamp groups
  # Use timestamp_group for selection labels; backed by dbscan_cluster
  # Ensure timestamp_group exists (created in main.R)
  if ('timestamp_group' %in% names(merged_table)) {
    timestamp_groups <- getUnique(merged_table, 'timestamp_group')
  } else if ('timestamp' %in% names(merged_table)) {
    # Derive groups directly from timestamps when group column is absent
    timestamp_groups <- sort(unique(merged_table$timestamp))
  } else {
    timestamp_groups <- as.POSIXct(character(0))
  }
  # Derive cluster label per timestamp_group to avoid length mismatch (only if both columns exist)
  if (all(c('timestamp_group','dbscan_cluster') %in% names(merged_table))) {
    tg_df <- merged_table %>%
      dplyr::select(timestamp_group, dbscan_cluster) %>%
      dplyr::distinct() %>%
      dplyr::arrange(timestamp_group)
    # Align clusters to timestamp_groups order
    cluster_labels <- tg_df %>%
      dplyr::group_by(timestamp_group) %>%
      dplyr::summarise(dbscan_cluster = dplyr::first(dbscan_cluster), .groups = 'drop') %>%
      dplyr::arrange(timestamp_group) %>%
      dplyr::pull(dbscan_cluster)
    named_timestamp_groups <- timestamp_groups
    names(named_timestamp_groups) <- paste(
      format(timestamp_groups, format = '%m-%d %H:%M'),
      sprintf('(cluster %s)', as.character(cluster_labels))
    )
  } else {
    # Fallback: labels without cluster info
    named_timestamp_groups <- timestamp_groups
    if (length(named_timestamp_groups) > 0) {
      names(named_timestamp_groups) <- format(timestamp_groups, format = '%m-%d %H:%M')
    }
  }
  
  # Named timestamp groups prepared above
  
  # Output variables
  numeric_cols <- names(merged_table)[sapply(merged_table, is.numeric)]
  out_variables <- sort(numeric_cols)
  
  return(list(
    anova_factors = anova_factors,
    treatments = treatments,
    cultivars = cultivars,
    timestamp_groups = timestamp_groups,
    named_timestamp_groups = named_timestamp_groups,
    out_variables = out_variables
  ))
}

# =============================================================================
# DRY UTILITIES: LOGIT AND DBSCAN
# =============================================================================

#' Apply logit transform to percentage columns and optionally treat infinities
#' @param data Data frame
#' @param use_logit logical: whether to apply transform
#' @param treat_inf logical: replace ±Inf with nearest finite bounds
#' @param selected_columns character: specific columns to transform (optional, defaults to all percentage columns)
#' @return Mutated data frame
applyLogitTransform <- function(data, use_logit, treat_inf = TRUE, selected_columns = NULL) {
  if (!isTRUE(use_logit)) return(data)
  fix_perc_imprecision <- function(x) {
    dplyr::case_when(
      (x >= 0) & (x <= 1.00) ~ x,
      (x < 0) & (x >= -0.01) ~ 0,
      (x > 1) & (x <= 1.01) ~ 1,
      .default = NA
    )
  }
  
  # Determine which columns to transform
  if (!is.null(selected_columns)) {
    # Transform only selected columns that are percentage columns
    cols_to_transform <- selected_columns[grepl('_percent$', selected_columns) & selected_columns %in% names(data)]
  } else {
    # Transform all percentage columns (original behavior)
    cols_to_transform <- names(data)[grepl('_percent$', names(data))]
  }
  
  if (length(cols_to_transform) == 0) return(data)
  
  # Apply transformation to selected columns
  for (col in cols_to_transform) {
    logit_col <- stringi::stri_replace_all_fixed(col, pattern = '_percent', replacement = '_logit')
    data[[logit_col]] <- stats::qlogis(fix_perc_imprecision(data[[col]]))
  }
  if (isTRUE(treat_inf)) {
    # Handle infinities for the transformed columns
    for (col in cols_to_transform) {
      logit_col <- stringi::stri_replace_all_fixed(col, pattern = '_percent', replacement = '_logit')
      finite_vals <- data[[logit_col]][!is.na(data[[logit_col]]) & !is.infinite(data[[logit_col]])]
      if (length(finite_vals) > 0) {
        min_val <- floor(min(finite_vals))
        max_val <- ceiling(max(finite_vals))
        data[[logit_col]] <- dplyr::case_when(
          is.infinite(data[[logit_col]]) & data[[logit_col]] < 0 ~ min_val,
          is.infinite(data[[logit_col]]) & data[[logit_col]] > 0 ~ max_val,
          TRUE ~ data[[logit_col]]
        )
      }
    }
  }
  data
}

#' Compute DBSCAN clusters from timestamps and eps in hours
#' @param timestamps POSIXct vector
#' @param eps_hours numeric eps in hours
#' @return integer cluster vector
computeDbscanClustersFromTimestamps <- function(timestamps, eps_hours) {
  # Validate eps parameter
  if (is.null(eps_hours) || is.na(eps_hours) || eps_hours < 0) {
    logError(paste('Invalid eps_hours parameter:', eps_hours))
    return(rep(1, length(timestamps)))  # Return single cluster as fallback
  }
  
  hours_from_start <- as.numeric(difftime(timestamps, min(timestamps), units = 'hours'))
  dbscan::dbscan(matrix(hours_from_start, ncol = 1), eps = eps_hours)$cluster
}

#' Preprocess data: fix single-level factors and ensure consistent data types
#' @param data data.frame to preprocess
#' @return preprocessed data.frame
preprocessDataForAnalysis <- function(data) {
  # Ensure dplyr is available
  if (!requireNamespace('dplyr', quietly = TRUE)) {
    stop('dplyr package is required for data preprocessing')
  }
  # Ensure critical metadata columns are always character type
  metadata_cols <- c('genotype', 'g_alias', 'cultivar', 'treatment')
  for (col in metadata_cols) {
    if (col %in% names(data)) {
      data[[col]] <- as.character(data[[col]])
    }
  }
  
  # Identify and drop single-level factors (except required ones)
  required_factors <- c('timestamp', 'unit', 'dbscan_cluster', 'timestamp_group')
  factor_cols <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]
  factor_cols <- factor_cols[!factor_cols %in% required_factors]
  
  single_level_factors <- c()
  for (col in factor_cols) {
    unique_vals <- length(unique(data[[col]][!is.na(data[[col]])]))
    if (unique_vals <= 1) {
      single_level_factors <- c(single_level_factors, col)
      cat("⚠️ Dropping single-level factor:", col, "(unique values:", unique_vals, ")\n")
    }
  }
  
  # Drop single-level factors
  if (length(single_level_factors) > 0) {
    data <- data %>%
      dplyr::select(-all_of(single_level_factors))
    cat("✅ Dropped", length(single_level_factors), "single-level factors:", paste(single_level_factors, collapse = ", "), "\n")
  }
  
  # Apply droplevels to all factor columns to clean up unused levels
  factor_cols_remaining <- names(data)[sapply(data, function(x) is.factor(x))]
  for (col in factor_cols_remaining) {
    data[[col]] <- droplevels(data[[col]])
  }
  
  return(data)
}

#' Ensure dbscan_cluster on a data frame given timestamp and eps
#' @param df Data frame with timestamp
#' @param eps_hours numeric eps hours
#' @return Data frame with factor dbscan_cluster
ensureDbscanOnFrame <- function(df, eps_hours) {
  if (!('timestamp' %in% names(df))) return(df)
  cl <- computeDbscanClustersFromTimestamps(df$timestamp, eps_hours)
  df$dbscan_cluster <- as.factor(cl)
  df
}

#' Renumber clusters by ascending timeline (earliest cluster -> 1, next -> 2, ...)
#' @param df Data frame with timestamp and dbscan_cluster
#' @return Data frame with dbscan_cluster remapped to ordered numeric factor
renumberClustersByTime <- function(df) {
  if (!all(c('timestamp','dbscan_cluster') %in% names(df))) return(df)
  # Compute earliest timestamp per cluster
  order_map <- df %>%
    dplyr::group_by(dbscan_cluster) %>%
    dplyr::summarise(tmin = min(timestamp, na.rm = TRUE), .groups = 'drop') %>%
    dplyr::arrange(tmin) %>%
    dplyr::mutate(new_id = dplyr::row_number()) %>%
    dplyr::select(dbscan_cluster, new_id)
  # Join back and relabel
  df <- df %>% dplyr::left_join(order_map, by = 'dbscan_cluster')
  df$dbscan_cluster <- factor(df$new_id, levels = sort(unique(df$new_id)), ordered = TRUE)
  df$new_id <- NULL
  df
}

