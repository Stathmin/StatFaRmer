# StatFaRmer Minimal Deployment

This is a minimal deployment version of StatFaRmer that uses the existing infrastructure but pre-selects the `project_NO3` project for ShinyApps.io deployment.

## What's Different

- **Pre-selected project**: Automatically loads `project_NO3` (wheat NO3 data)
- **No project selection UI**: Skips the project selection step
- **Full functionality**: Uses all existing modules and features
- **Minimal changes**: Reuses existing code with minimal modifications

## Files Structure

```
deploy_minimal/
├── app.R              # Main application file
├── deploy.R           # Deployment script with package installation
├── README.md          # This file
├── src/               # Symlink to ../src
├── shiny/             # Symlink to ../shiny
├── config/            # Symlink to ../config
└── data/              # Symlink to ../data
```

## Deployment to ShinyApps.io

### Option 1: Using rsconnect

```r
# Install rsconnect if needed
install.packages("rsconnect")

# Set up account (one time)
rsconnect::setAccountInfo(
  name = "your-account-name",
  token = "your-token",
  secret = "your-secret"
)

# Deploy
rsconnect::deployApp(
  appDir = "deploy_minimal",
  appName = "statfarmer-minimal",
  appFile = "app.R",
  forceUpdate = TRUE
)
```

### Option 2: Manual Upload

1. Zip the entire `deploy_minimal` folder
2. Upload to ShinyApps.io via web interface
3. Use `app.R` as the main file

## Features

- **Full StatFaRmer functionality** with all modules
- **Pre-loaded project_NO3 data** (wheat NO3 analysis)
- **All visualization options** (time series, box plots, etc.)
- **Complete statistical analysis** (ANOVA, mixed-effects, post-hoc)
- **Effect sizes and growth summaries**
- **Data and plot export**
- **No project selection needed** - ready to use immediately

## Data

The app automatically loads the wheat NO3 project data including:
- 3,915 observations across 58 variables
- 45+ phenotypic variables (biomass, height, leaf area, etc.)
- Treatment and cultivar factors
- Pre-calculated outlier detection and clustering
- Multiple timestamp groups for time series analysis

## Advantages

- ✅ **Minimal code changes** - reuses existing infrastructure
- ✅ **Full functionality** - no feature concessions
- ✅ **Easy deployment** - just pre-selects a project
- ✅ **Maintainable** - uses existing modules and utilities
- ✅ **No data duplication** - uses symlinks to existing data

## Customization

To use with a different project:
1. Change `options(statfarmer.project = "project_NAME")` in `app.R`
2. Ensure the project exists in the `data/` directory
3. Redeploy

## Performance

- **Startup time**: ~15-20 seconds (full app loading)
- **Memory usage**: ~300-500 MB
- **Response time**: <3 seconds for most operations
- **File size**: ~100-200 MB total (including all dependencies)
