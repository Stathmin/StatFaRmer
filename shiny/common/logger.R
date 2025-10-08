# Simple logging helpers for Shiny

if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)

ensureLogDirs <- function() {
  logs_dir <- here('logs')
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
}

currentLogLevel <- function() {
  lvl <- Sys.getenv('LOG_LEVEL', unset = 'INFO')
  toupper(lvl)
}

levels_order <- c('DEBUG' = 10, 'INFO' = 20, 'WARN' = 30, 'ERROR' = 40)

shouldLog <- function(level) {
  levels_order[[toupper(level)]] >= levels_order[[currentLogLevel()]]
}

logInfo <- function(message) {
  ensureLogDirs()
  if (shouldLog('INFO')) {
    entry <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 'INFO', message, sep = ' | ')
    write(entry, file = here('logs', 'app.log'), append = TRUE)
  }
}

logWarn <- function(message) {
  ensureLogDirs()
  if (shouldLog('WARN')) {
    entry <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 'WARN', message, sep = ' | ')
    write(entry, file = here('logs', 'app.log'), append = TRUE)
  }
}

logError <- function(message) {
  ensureLogDirs()
  entry <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 'ERROR', message, sep = ' | ')
  write(entry, file = here('logs', 'error.log'), append = TRUE)
}

logDebug <- function(message) {
  ensureLogDirs()
  if (shouldLog('DEBUG')) {
    entry <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 'DEBUG', message, sep = ' | ')
    write(entry, file = here('logs', 'app.log'), append = TRUE)
  }
}

safe_to_json <- function(x) {
  if (!requireNamespace('jsonlite', quietly = TRUE)) return('')
  tryCatch(jsonlite::toJSON(x, auto_unbox = TRUE, null = 'null'), error = function(e) '')
}

logEvent <- function(level = 'INFO', event = 'event', data = list()) {
  ensureLogDirs()
  if (!shouldLog(level)) return(invisible(FALSE))
  payload <- safe_to_json(list(event = event, data = data))
  entry <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), toupper(level), payload, sep = ' | ')
  write(entry, file = here('logs', 'app.log'), append = TRUE)
  
  # Also write ERROR level events to error.log
  if (toupper(level) == 'ERROR') {
    write(entry, file = here('logs', 'error.log'), append = TRUE)
  }
  
  invisible(TRUE)
}

withBenchmark <- function(operation_name, expr) {
  bm <- startBenchmark(operation_name)
  on.exit({
    endBenchmark(bm)
  }, add = TRUE)
  force(expr)
}


