#!/usr/bin/env Rscript
# Headless main-app launcher for puppeteer smoke (does not edit launch_statfarmer.R).
# Usage:
#   Rscript scripts/browser/run_app_smoke.R [project_id] [port]
suppressPackageStartupMessages({
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  library(here)
  library(shiny)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) args[[1L]] else "project_2025-07-14_b51-91_soy"
port <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3848L

Sys.setenv(SELECTED_PROJECT = project)
options(statfarmer.project = project)
options(shiny.port = port)
options(shiny.host = "127.0.0.1")

cat(sprintf("Smoke app: project=%s port=%s\n", project, port))
runApp(here("shiny", "app"), port = port, host = "127.0.0.1", launch.browser = FALSE)
