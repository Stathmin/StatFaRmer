# StatFaRmer testthat entrypoint (project is not an R package).
# Run from repo root: Rscript tests/testthat.R
# Or a single file: Rscript -e "testthat::test_file('tests/testthat/test-fixture-mini.R')"

library(testthat)

test_dir(
  path = file.path("tests", "testthat"),
  reporter = "progress",
  stop_on_failure = FALSE
)
