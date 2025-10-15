# StatFaRmer Main Application User Guide

## Overview
The StatFaRmer Main Application provides statistical analysis, interactive visualizations, and export capabilities for wizard-processed plant phenotyping data.

## Launching the Application

```bash
# Launch main app directly
Rscript launch_statfarmer.R app

# Launch with specific project
Rscript launch_statfarmer.R app project_NO3

# Launch from wizard (recommended)
Rscript launch_statfarmer.R wizard
# Then click "Launch Main App" after processing
```

## Interface
- **Left Sidebar**: Project selection, variable selection, analysis parameters
- **Main Panel**: Results display with tabbed interface
- **Professional Theme**: cowplot themes for publication-ready visualizations

## Statistical Analysis Features

### 1. ANOVA Analysis

**Automatic Model Selection**
- **Classic ANOVA**: For balanced designs with ≤10 factor levels
- **Mixed Models**: For unbalanced designs or 11-20 factor levels  
- **Spline Models**: For time-series analysis with temporal effects

**Advanced Model Options**
- **Auto Selection**: Intelligent model choice based on factor cardinality
- **Force Options**: Override automatic selection (ANOVA/LMM/Spline)
- **Formula Preview**: Real-time preview with cardinality warnings
- **Model Info Box**: Diagnostics, warnings, and model details

**Model Features**
- **Factor Interactions**: Automatic two-way interactions
- **Random Effects**: Automatic detection of blocking factors (`dbscan_cluster`, `timestamp_group`)
- **Assumption Testing**: Normality, homogeneity, independence checks
- **Diagnostic Plots**: Residual analysis and QQ plots

### 2. Mixed Effects Models

**Automatic Detection**
- **High Cardinality Factors**: >10 levels trigger mixed model
- **Blocking Factors**: `dbscan_cluster`, `timestamp_group` as random effects
- **Performance Optimization**: Automatic fallback for large datasets

**Model Types**
- **lme4 Integration**: Full mixed-effects modeling
- **Random Intercepts**: `(1|blocking_factor)` structure
- **Fixed Effects**: Treatment and cultivar effects
- **Time Effects**: Temporal clustering as random effect

### 3. Spline Analysis

**Time-Series Modeling**
- **Temporal Clustering**: DBSCAN clusters as time points
- **Smooth Trends**: Spline fitting for growth curves
- **Treatment Comparison**: Compare growth patterns across treatments
- **Peak Analysis**: Identify maximum growth periods

## Visualization System

**Interactive Plots**
- **Plotly Integration**: Zoom, pan, hover information, legend controls
- **Plot Types**: ANOVA boxplots, time series, residual plots, QQ plots
- **Export Options**: High-resolution downloads (SVG, PNG, PDF)

**Professional Styling**
- **cowplot Themes**: Publication-ready appearance
- **Customization**: Facet formulas, color schemes, custom sizing
- **Resolution**: 300 DPI for publications

## Data Management

**Project Selection**
- **Wizard Integration**: Automatic project detection and configuration preservation
- **Data Validation**: Ensures processed data integrity
- **Project Types**: Standard projects, restricted projects (public deployment)

**Data Filtering**
- **Subsetting Options**: Treatment selection, cultivar selection, time clusters, outlier handling
- **Dynamic Updates**: Real-time filter application
- **Combined Filters**: Multiple criteria simultaneously

## Results & Export

**Analysis Results**
- **ANOVA Tables**: Detailed statistical results with model information
- **Descriptive Statistics**: Summary statistics by group
- **Post-hoc Tests**: Enhanced Tukey HSD with robust letters generation
- **Group Letters**: Compact letter display with intelligent filtering
- **Model Diagnostics**: Assumption testing results and convergence checks

**Enhanced Post-hoc Analysis**
- **Robust Letters Generation**: Multiple fallback strategies for CLD generation
- **Intelligent Filtering**: Automatic filtering for large comparison sets (1000+)
- **Stratified Analysis**: Support for by/stratum columns in letters
- **Error Recovery**: Graceful degradation when analysis fails

**Export Capabilities**
- **Plot Downloads**: High resolution (300 DPI), multiple formats, custom sizing
- **Data Exports**: Filtered data, results tables, summary reports, configuration
- **Letters Tables**: Downloadable group letters with significance information

## Specialized Modules

### Effect Sizes Module
- **Cohen's d Calculations**: Standardized differences with confidence intervals
- **Practical Significance**: Beyond statistical significance
- **Comparison Tables**: Effect sizes across treatments

### Growth Summaries Module
- **Temporal Analysis**: Area Under Curve (AUC), peak analysis, growth rates
- **Timing Metrics**: Developmental milestones

## Troubleshooting

**Data Loading Problems**
- **Missing Files**: Ensure wizard processing completed
- **File Permissions**: Check data directory access
- **Configuration Errors**: Verify config.json format

**Analysis Errors**
- **Model Failures**: Check factor cardinality and data structure
- **Convergence Issues**: Adjust model parameters or data
- **Memory Limits**: Use median aggregation for large datasets

**Performance Optimization**
- **Large Datasets**: Use median aggregation, reduce factor levels, filter data
- **Model Selection**: Start simple, progressive complexity, check assumptions

---

**Data Preparation**: Use the [Master Wizard](wizard_guide.md) to prepare your data before analysis.

**Statistical Methods**: See [Statistical Methods Reference](common_guide.md) for technical details on model selection and analysis approaches.