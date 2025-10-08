# StatFaRmer Data Validation Module
# Iteration 2: Data validation

library(tidyverse)
library(readxl)
library(zip)

# =============================================================================
# MAIN VALIDATION FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# Schema normalization helpers
# -----------------------------------------------------------------------------

#' Normalize metadata column names to expected schema
#' @param data Data frame
#' @return Normalized data frame
normalizeMetadataSchema <- function(data) {
  colmap <- c(
    'v.t.r' = 'V.T.R',
    'vtr' = 'V.T.R',
    'plant' = 'Plant',
    'var' = 'Plant',
    'treatment' = 'Treatment',
    'repeat' = 'Repeat',
    'cultivar' = 'Cultivar'
  )
  current <- colnames(data)
  lower <- tolower(current)
  new_names <- current
  for (i in seq_along(lower)) {
    if (lower[i] %in% names(colmap)) new_names[i] <- colmap[[lower[i]]]
  }
  colnames(data) <- new_names
  # type normalization
  if ("Plant" %in% names(data)) data$Plant <- suppressWarnings(as.numeric(data$Plant))
  if ("Treatment" %in% names(data)) data$Treatment <- suppressWarnings(as.numeric(data$Treatment))
  if ("Repeat" %in% names(data)) data$Repeat <- suppressWarnings(as.numeric(data$Repeat))
  if ("Cultivar" %in% names(data)) data$Cultivar <- as.character(data$Cultivar)
  return(data)
}

#' Normalize translation column names to expected schema
#' @param data Data frame
#' @return Normalized data frame
normalizeTranslationSchema <- function(data) {
  colmap <- c(
    't:x:y' = 'T:X:Y',
    't.x.y' = 'T:X:Y',
    'v.t.r' = 'V.T.R',
    'vtr' = 'V.T.R'
  )
  current <- colnames(data)
  lower <- tolower(current)
  new_names <- current
  for (i in seq_along(lower)) {
    if (lower[i] %in% names(colmap)) new_names[i] <- colmap[[lower[i]]]
  }
  colnames(data) <- new_names
  return(data)
}

#' Validate TraitFinder data (ZIP files)
#' @param zip_path Path to ZIP file with TraitFinder data
#' @return List with validation results
validateTraitFinderData <- function(zip_path) {
  result <- list(
    valid = FALSE,
    errors = character(0),
    warnings = character(0),
    info = character(0)
  )
  
  tryCatch({
    # Check file existence
    if (!file.exists(zip_path)) {
      result$errors <- c(result$errors, paste("File not found:", zip_path))
      return(result)
    }
    
    # Check file extension
    if (!str_detect(zip_path, "\\.zip$")) {
      result$errors <- c(result$errors, "File must have .zip extension")
      return(result)
    }
    
    # Check ZIP contents
    zip_contents <- zip_list(zip_path)
    if (nrow(zip_contents) == 0) {
      result$errors <- c(result$errors, "ZIP file is empty")
      return(result)
    }
    
    # Check for CSV files
    csv_files <- zip_contents$filename[str_detect(zip_contents$filename, "\\.csv$")]
    if (length(csv_files) == 0) {
      result$errors <- c(result$errors, "No CSV files found in ZIP")
      return(result)
    }
    
    result$info <- c(result$info, paste("Found CSV files:", length(csv_files)))
    result$valid <- TRUE
    
  }, error = function(e) {
    result$errors <- c(result$errors, paste("Error reading ZIP:", e$message))
  })
  
  return(result)
}

#' Validate metadata (handmade CSV)
#' @param csv_path Path to CSV file with metadata
#' @return List with validation results
validateMetadata <- function(csv_path) {
  result <- list(
    valid = FALSE,
    errors = character(0),
    warnings = character(0),
    info = character(0)
  )
  
  tryCatch({
    # Check file existence
    if (!file.exists(csv_path)) {
      result$errors <- c(result$errors, paste("File not found:", csv_path))
      return(result)
    }
    
    # Read file
    data <- read_csv(csv_path, show_col_types = FALSE)
    data <- normalizeMetadataSchema(data)
    
    # Required fields
    required_fields <- c("V.T.R", "Plant", "Treatment", "Repeat", "Cultivar")
    missing_fields <- setdiff(required_fields, colnames(data))
    
    if (length(missing_fields) > 0) {
      result$errors <- c(result$errors, 
                        paste("Missing required fields:", 
                              paste(missing_fields, collapse = ", ")))
    }
    
    # Check data types
    if ("Plant" %in% colnames(data)) {
      if (!is.numeric(data$Plant)) {
        result$warnings <- c(result$warnings, "Plant field should be numeric")
      }
    }
    
    if ("Treatment" %in% colnames(data)) {
      if (!is.numeric(data$Treatment)) {
        result$warnings <- c(result$warnings, "Treatment field should be numeric")
      }
    }
    
    if ("Repeat" %in% colnames(data)) {
      if (!is.numeric(data$Repeat)) {
        result$warnings <- c(result$warnings, "Repeat field should be numeric")
      }
    }
    
    # Check V.T.R uniqueness
    if ("V.T.R" %in% colnames(data)) {
      duplicates <- sum(duplicated(data$V.T.R))
      if (duplicates > 0) {
        result$warnings <- c(result$warnings, 
                            paste("Found V.T.R duplicates:", duplicates))
      }
    }
    
    # Check for empty values
    empty_rows <- sum(apply(data, 1, function(x) any(is.na(x) | x == "")))
    if (empty_rows > 0) {
      result$warnings <- c(result$warnings, 
                          paste("Found rows with empty values:", empty_rows))
    }
    
    result$info <- c(result$info, 
                    paste("Data rows:", nrow(data)),
                    paste("Columns:", ncol(data)))
    
    if (length(result$errors) == 0) {
      result$valid <- TRUE
    }
    
  }, error = function(e) {
    result$errors <- c(result$errors, paste("Error reading CSV:", e$message))
  })
  
  return(result)
}

#' Validate coordinates (translation CSV)
#' @param csv_path Path to CSV file with coordinates
#' @return List with validation results
validateCoordinates <- function(csv_path) {
  result <- list(
    valid = FALSE,
    errors = character(0),
    warnings = character(0),
    info = character(0)
  )
  
  tryCatch({
    # Check file existence
    if (!file.exists(csv_path)) {
      result$errors <- c(result$errors, paste("File not found:", csv_path))
      return(result)
    }
    
    # Read file
    data <- read_csv(csv_path, show_col_types = FALSE)
    data <- normalizeTranslationSchema(data)
    
    # Required fields
    required_fields <- c("T:X:Y", "V.T.R")
    missing_fields <- setdiff(required_fields, colnames(data))
    
    if (length(missing_fields) > 0) {
      result$errors <- c(result$errors, 
                        paste("Missing required fields:", 
                              paste(missing_fields, collapse = ", ")))
    }
    
    # Check T:X:Y coordinate format
    if ("T:X:Y" %in% colnames(data)) {
      coord_pattern <- "^\\d+:\\d+:\\d+$"
      invalid_coords <- !str_detect(data$`T:X:Y`, coord_pattern)
      
      if (any(invalid_coords)) {
        invalid_count <- sum(invalid_coords)
        result$warnings <- c(result$warnings, 
                            paste("Found invalid coordinates:", invalid_count))
      }
    }
    
    # Check V.T.R uniqueness
    if ("V.T.R" %in% colnames(data)) {
      duplicates <- sum(duplicated(data$V.T.R))
      if (duplicates > 0) {
        result$warnings <- c(result$warnings, 
                            paste("Found V.T.R duplicates:", duplicates))
      }
    }
    
    # Check for empty values
    empty_rows <- sum(apply(data, 1, function(x) any(is.na(x) | x == "")))
    if (empty_rows > 0) {
      result$warnings <- c(result$warnings, 
                          paste("Found rows with empty values:", empty_rows))
    }
    
    result$info <- c(result$info, 
                    paste("Data rows:", nrow(data)),
                    paste("Columns:", ncol(data)))
    
    if (length(result$errors) == 0) {
      result$valid <- TRUE
    }
    
  }, error = function(e) {
    result$errors <- c(result$errors, paste("Error reading CSV:", e$message))
  })
  
  return(result)
}

#' Validate groups (Excel file)
#' @param excel_path Path to Excel file with groups
#' @return List with validation results
validateGroups <- function(excel_path) {
  result <- list(
    valid = FALSE,
    errors = character(0),
    warnings = character(0),
    info = character(0)
  )
  
  tryCatch({
    # Check file existence
    if (!file.exists(excel_path)) {
      result$errors <- c(result$errors, paste("File not found:", excel_path))
      return(result)
    }
    
    # Check file extension
    if (!str_detect(excel_path, "\\.(xlsx|xls)$")) {
      result$errors <- c(result$errors, "File must have .xlsx or .xls extension")
      return(result)
    }
    
    # Read file
    data <- read_excel(excel_path)
    
    # Check for empty data
    if (nrow(data) == 0) {
      result$errors <- c(result$errors, "Excel file is empty")
      return(result)
    }
    
    # Check for empty values
    empty_rows <- sum(apply(data, 1, function(x) any(is.na(x) | x == "")))
    if (empty_rows > 0) {
      result$warnings <- c(result$warnings, 
                          paste("Found rows with empty values:", empty_rows))
    }
    
    result$info <- c(result$info, 
                    paste("Data rows:", nrow(data)),
                    paste("Columns:", ncol(data)))
    
    result$valid <- TRUE
    
  }, error = function(e) {
    result$errors <- c(result$errors, paste("Error reading Excel:", e$message))
  })
  
  return(result)
}

# =============================================================================
# INTEGRATED PROJECT VALIDATION
# =============================================================================

#' Complete project validation
#' @param project_path Path to project folder
#' @return List with validation results for all files
validateProject <- function(project_path) {
  result <- list(
    valid = FALSE,
    errors = character(0),
    warnings = character(0),
    info = character(0),
    files = list()
  )
  
  # Find files in project
  zip_files <- list.files(project_path, pattern = "*_data\\.zip$", full.names = TRUE)
  handmade_files <- list.files(project_path, pattern = "*_handmade\\.csv$", full.names = TRUE)
  translation_files <- list.files(project_path, pattern = "*_translation\\.csv$", full.names = TRUE)
  groups_files <- list.files(project_path, pattern = "groups\\.xlsx$", full.names = TRUE)
  
  # Validate ZIP files
  for (zip_file in zip_files) {
    result$files[[basename(zip_file)]] <- validateTraitFinderData(zip_file)
  }
  
  # Validate handmade files
  for (handmade_file in handmade_files) {
    result$files[[basename(handmade_file)]] <- validateMetadata(handmade_file)
  }
  
  # Validate translation files
  for (translation_file in translation_files) {
    result$files[[basename(translation_file)]] <- validateCoordinates(translation_file)
  }
  
  # Validate groups files
  for (groups_file in groups_files) {
    result$files[[basename(groups_file)]] <- validateGroups(groups_file)
  }
  
  # Aggregate results
  for (file_result in result$files) {
    result$errors <- c(result$errors, file_result$errors)
    result$warnings <- c(result$warnings, file_result$warnings)
    result$info <- c(result$info, file_result$info)
  }
  
  # Overall validation status
  all_valid <- all(sapply(result$files, function(x) x$valid))
  result$valid <- all_valid && length(result$errors) == 0
  
  return(result)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Format validation results for display
#' @param validation_result Validation result
#' @return Formatted string
formatValidationResult <- function(validation_result) {
  if (validation_result$valid) {
    status <- "✅ VALIDATION PASSED"
  } else {
    status <- "❌ VALIDATION FAILED"
  }
  
  output <- c(status, "")
  
  if (length(validation_result$errors) > 0) {
    output <- c(output, "🚨 ERRORS:", validation_result$errors, "")
  }
  
  if (length(validation_result$warnings) > 0) {
    output <- c(output, "⚠️ WARNINGS:", validation_result$warnings, "")
  }
  
  if (length(validation_result$info) > 0) {
    output <- c(output, "ℹ️ INFO:", validation_result$info, "")
  }
  
  return(paste(output, collapse = "\n"))
}

#' Log validation results
#' @param validation_result Validation result
#' @param project_name Project name
logValidationResult <- function(validation_result, project_name) {
  # Create logs directory if it doesn't exist
  if (!dir.exists("logs")) {
    dir.create("logs", recursive = TRUE)
  }
  
  # Write to log
  log_entry <- paste(
    Sys.time(),
    "VALIDATION",
    project_name,
    ifelse(validation_result$valid, "SUCCESS", "FAILED"),
    paste(validation_result$errors, collapse = "; "),
    sep = " | "
  )
  
  write(log_entry, file = "logs/validation.log", append = TRUE)
}
