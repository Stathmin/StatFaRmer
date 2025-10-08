# Function to check system dependencies
checkSystemDeps <- function() {
  message("Checking system dependencies...")
  
  # Check if we're on Linux
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Linux") {
    message("Linux detected. Some packages may need additional system libraries.")
    message("If installation fails, you may need to install:")
    message("  sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libcairo2-dev libharfbuzz-dev libfribidi-dev")
  }
}

# Function to check if renv is installed and install it if not
checkRenv <- function() {
  if (!requireNamespace('renv', quietly = TRUE)) {
    message("Installing renv package...")
    # Install to user library (R's standard approach)
    install.packages('renv', repos = 'https://cran.rstudio.com/')
    message("renv installed successfully!")
  } else {
    message("renv is already installed.")
  }
}

# Function to ensure renv environment is properly set up
setupRenv <- function() {
  # Set package installation preferences based on platform
  if (.Platform$OS.type == 'windows') {
    # Windows: prefer binary packages for speed
    options(install.packages.compile.from.source = "ifneeded")
    options(install.packages.check.source = "yes")
    install_type <- "binary"
  } else {
    # Linux: use source packages
    install_type <- "source"
  }
  
  # Check if renv is already initialized
  if (!file.exists('renv.lock')) {
    message("renv.lock not found. Initializing renv environment...")
    renv::init()
  } else {
    message("renv.lock found. Restoring environment...")
    # Restore packages (renv handles platform-specific installation automatically)
    renv::restore(prompt = FALSE)
  }
}

# Function to install failed packages with better error handling
installFailedPackages <- function() {
  message("Attempting to install previously failed packages...")
  
  # List of packages that commonly fail due to system dependencies
  failed_packages <- c(
    "curl", "httr", "systemfonts", "textshaping", "ragg",
    "svglite", "gdtools", "officer", "flextable", "tidyverse",
    "car"
  )
  
  for (pkg in failed_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      tryCatch({
        install.packages(pkg, repos = 'https://cran.rstudio.com/')
        message(paste(pkg, "installed successfully!"))
      }, error = function(e) {
        message(paste("Failed to install", pkg, ":", e$message))
        message("This package may need additional system libraries.")
      })
    }
  }
}

installStatfarmer <- function() {
  message("Starting installation of StatFaRmer...")

  # Check system dependencies first
  checkSystemDeps()
  
  message("Checking for renv package...")
  checkRenv()

  message("Setting up renv environment...")
  setupRenv()
  
  # Try to install failed packages
  installFailedPackages()

  message("Installation of StatFaRmer completed!")
  message("You can now run: Rscript src/main.R")
}

installStatfarmer()
