# Guard: tidyverse attaches jsonlite, which masks shiny::validate.
# Bare validate(need(...), need(...)) then fails with "unused arguments".

test_that("BLUP/DR modules use shiny::validate (not bare validate)", {
  root <- here::here()
  files <- c(
    file.path(root, "shiny/app/module_blup_ranking.R"),
    file.path(root, "shiny/app/module_trait_reduction.R")
  )
  for (f in files) {
    expect_true(file.exists(f), info = f)
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    bare <- grepl("(^|\n)[[:space:]]*validate\\(", txt, perl = TRUE)
    expect_false(bare, info = paste("bare validate( in", basename(f)))
    expect_true(grepl("shiny::validate\\(", txt), info = paste("missing shiny::validate in", basename(f)))
  }
})

test_that("jsonlite::validate rejects multi-need args (user crash)", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("shiny")
  validate <- jsonlite::validate
  err <- tryCatch(
    validate(shiny::need(TRUE, "a"), shiny::need(TRUE, "b")),
    error = function(e) e
  )
  expect_true(inherits(err, "error"))
  expect_match(conditionMessage(err), "unused argument|неиспользован", ignore.case = TRUE)
  expect_error(
    shiny::validate(shiny::need(FALSE, "blocked")),
    class = "shiny.silent.error"
  )
})
