# Function to check if renv is installed and install it if not
check_renv <- function() {
  if (!requireNamespace('renv', quietly = TRUE)) {
    install.packages('renv')
  }
}

# Function to install flextable package
install_flextable <- function() {
  install.packages('flextable',
                   type = if (.Platform$OS.type == 'windows') 'binary' else 'source')
}

# Main installation function
install_statfarmer <- function() {
  # Print starting message
  message("Starting installation of StatFaRmer...")

  # Check and install renv if necessary
  message("Checking for renv package...")
  check_renv()

  # Initialize and hydrate renv environment
  message("Initializing renv environment...")
  renv::init()

  message("Hydrating renv environment...")
  renv::install(prompt = FALSE)

  # Install flextable package
  message("Installing flextable package...")
  install_flextable()

  # Print completion message
  message("Installation of StatFaRmer completed successfully!")
}

# Run the installation function
install_statfarmer()
