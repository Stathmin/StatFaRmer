# StatFaRmer Minimal Deployment Launcher
# Simple launcher for local testing

# Set working directory to project root
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Source the minimal deployment app
source('deploy_minimal/app.R')
