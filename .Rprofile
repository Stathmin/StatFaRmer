# StatFaRmer R Profile
# Ensure proper library paths for Cursor R Tools

# Set up renv if not already done
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# Ensure languageserver is available
if (!requireNamespace("languageserver", quietly = TRUE)) {
  install.packages("languageserver")
}

# Set options for better Cursor integration
options(
  repos = c(CRAN = "https://cloud.r-project.org/"),
  renv.config.auto.snapshot = TRUE,
  renv.config.ignored.packages = c("logger", "pracma", "vecsets", "writexl")
)