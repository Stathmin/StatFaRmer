# Shared test helpers for StatFaRmer (non-package project).
# Sourced automatically by testthat before tests.

sf_repo_root <- function() {
  if (requireNamespace("here", quietly = TRUE)) {
    return(here::here())
  }
  normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
}

#' Source a file under shiny/common/ without attaching Shiny.
#' @param name File name relative to shiny/common (e.g. "txy_parse.R")
#' @param envir Environment to source into (default: calling frame)
sf_source_common <- function(name, envir = parent.frame()) {
  path <- file.path(sf_repo_root(), "shiny", "common", name)
  if (!file.exists(path)) {
    stop("Missing shiny/common file: ", path, call. = FALSE)
  }
  sys.source(path, envir = envir)
}

#' Path to a committed test fixture under tests/testthat/fixtures/
sf_fixture_path <- function(...) {
  testthat::test_path("fixtures", ...)
}
