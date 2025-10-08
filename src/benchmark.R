# StatFaRmer Benchmarking Functions
# Iteration 3.5: Performance monitoring

# Load required libraries
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)

#' Start benchmark timer
#' @param operation_name Name of the operation being measured
#' @return Benchmark object
startBenchmark <- function(operation_name) {
  list(
    operation = operation_name,
    start_time = Sys.time(),
    start_memory = gc(verbose = FALSE)[2, 2] # Memory in MB
  )
}

#' End benchmark timer and log results
#' @param benchmark Benchmark object from startBenchmark
#' @param additional_info Additional information to log
#' @return Benchmark results
endBenchmark <- function(benchmark, additional_info = "") {
  end_time <- Sys.time()
  end_memory <- gc(verbose = FALSE)[2, 2]
  
  duration <- as.numeric(difftime(end_time, benchmark$start_time, units = "secs"))
  memory_used <- end_memory - benchmark$start_memory
  
  results <- list(
    operation = benchmark$operation,
    duration_seconds = round(duration, 3),
    memory_mb = round(memory_used, 2),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    additional_info = additional_info
  )
  
  # Log to benchmark file
  logBenchmark(results)
  
  return(results)
}

#' Log benchmark results to CSV file
#' @param results Benchmark results
logBenchmark <- function(results) {
  # Ensure logs directory exists
  logs_dir <- here('logs')
  if (!dir.exists(logs_dir)) {
    dir.create(logs_dir, recursive = TRUE)
  }
  
  benchmark_file <- here('logs', 'benchmark.csv')
  
  # Create data frame
  benchmark_df <- data.frame(
    timestamp = results$timestamp,
    operation = results$operation,
    duration_seconds = results$duration_seconds,
    memory_mb = results$memory_mb,
    additional_info = results$additional_info,
    stringsAsFactors = FALSE
  )
  
  # Append to file (create if doesn't exist)
  if (file.exists(benchmark_file)) {
    write.table(benchmark_df, file = benchmark_file, 
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  } else {
    write.csv(benchmark_df, file = benchmark_file, row.names = FALSE)
  }
  
  # Also print to console for immediate feedback
  cat(sprintf("⏱️  %s: %.3fs, %.2fMB - %s\n", 
              results$operation, results$duration_seconds, results$memory_mb, results$additional_info))
}

#' Get benchmark summary
#' @param operation_name Optional operation name to filter
#' @return Summary statistics
getBenchmarkSummary <- function(operation_name = NULL) {
  benchmark_file <- here('logs', 'benchmark.csv')
  
  if (!file.exists(benchmark_file)) {
    return(data.frame(message = "No benchmark data available"))
  }
  
  data <- read.csv(benchmark_file, stringsAsFactors = FALSE)
  
  if (!is.null(operation_name)) {
    data <- data[data$operation == operation_name, ]
  }
  
  if (nrow(data) == 0) {
    return(data.frame(message = "No data for specified operation"))
  }
  
  summary_stats <- data.frame(
    operation = if (is.null(operation_name)) "All Operations" else operation_name,
    total_operations = nrow(data),
    avg_duration = round(mean(data$duration_seconds), 3),
    min_duration = round(min(data$duration_seconds), 3),
    max_duration = round(max(data$duration_seconds), 3),
    total_duration = round(sum(data$duration_seconds), 3),
    avg_memory = round(mean(data$memory_mb), 2),
    max_memory = round(max(data$memory_mb), 2),
    stringsAsFactors = FALSE
  )
  
  return(summary_stats)
}

#' Check if performance meets targets
#' @param operation_name Operation to check
#' @param target_seconds Target duration in seconds
#' @return Performance check results
checkPerformanceTargets <- function(operation_name, target_seconds = 30) {
  summary_stats <- getBenchmarkSummary(operation_name)
  
  if ("message" %in% names(summary_stats)) {
    return(list(
      meets_target = FALSE,
      message = summary_stats$message
    ))
  }
  
  meets_target <- summary_stats$avg_duration <= target_seconds
  
  return(list(
    meets_target = meets_target,
    avg_duration = summary_stats$avg_duration,
    target_duration = target_seconds,
    performance_ratio = round(summary_stats$avg_duration / target_seconds, 2),
    message = if (meets_target) {
      sprintf("✅ Performance target met: %.3fs <= %.1fs", 
              summary_stats$avg_duration, target_seconds)
    } else {
      sprintf("⚠️  Performance target missed: %.3fs > %.1fs", 
              summary_stats$avg_duration, target_seconds)
    }
  ))
}







