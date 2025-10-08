# StatFaRmer Master Wizard User Guide

## Overview
The StatFaRmer Master Wizard prepares plant phenotyping data for statistical analysis through validation, preprocessing, and configuration management.

## Launching the Wizard

```bash
Rscript launch_statfarmer.R
# Select option [1] Launch Master Wizard
```

## Project Requirements

```
data/project_NAME/
├── *_data.zip              # TraitFinder experiment archive
├── *_handmade.csv          # Metadata with V.T.R, Treatment, Cultivar columns
├── *_translation.csv       # Spatial coordinates with V.T.R and T:X:Y columns
├── groups.xlsx             # Optional grouping factors
└── config.json             # Auto-generated configuration
```

## Workflow

### 1. Project Selection and Validation
- Select project from dropdown under "Input folder"
- Click "Validate input" to check file structure and data integrity
- Review validation results in the "Validation" tab

### 2. Processing Parameters

**DBSCAN Clustering**
- **DBSCAN eps (hours)**: Time window for clustering timestamps (default: 1 hour)
  - Smaller values: higher temporal resolution, more clusters
  - Larger values: lower temporal resolution, fewer clusters

**Data Transformations**
- **Apply logit transform**: Converts percentage data to logit scale (recommended for percentages)
- **Treat ±Inf/NaN**: Handles infinite values from logit transformation

**Data Preprocessing**
- **Single-level factor removal**: Automatically drops factors with only one unique value (prevents model fitting errors)
- **Metadata type consistency**: Ensures `genotype`, `g_alias`, `cultivar`, `treatment` are always character type
- **Factor level cleanup**: Removes unused factor levels with `droplevels()`

**Outlier Handling**
- **Enable outlier handling by cells**: Activates outlier detection
- **Outlier handling method**:
  - **Remove outliers**: Replaces outliers with NA (recommended)
  - **Winsorize outliers**: Caps outliers at threshold values
- **Detection method**:
  - **IQR**: Interquartile range method (robust, default)
  - **Z-score**: Standard deviation method (sensitive to distribution)

**Factor Configuration**
- **Factor columns for cells**: Select grouping factors for outlier detection
  - Outliers detected within each factor combination
  - Common factors: Treatment, Cultivar, Replication

**Technical Aggregation**
- **Aggregation of technical replicates**:
  - **Median**: Robust to outliers (recommended)
  - **Mean**: Standard average (sensitive to outliers)

**DBSCAN Outlier Clusters**
- **Select DBSCAN clusters as outliers**: Remove entire time clusters
  - Use for removing faulty measurement sessions
  - Preview clusters in "DBSCAN preview" tab

### 3. Preview and Processing
- Review parameters in "Preview params" tab
- Visualize clustering in "DBSCAN preview" tab
- Monitor logs in "Log" tab
- Click "Run processing and export" to execute pipeline

### 4. Launch Main Application
- After successful processing, click "Launch StatFaRmer"
- Or launch manually: `Rscript launch_statfarmer.R app project_NAME`

## Performance Benchmarks

**project_NO3 (58,380 rows, 49 columns)**
- Data Loading: 0.3-0.6 seconds (~25 MB memory)
- DBSCAN Clustering: 0.4-0.8 seconds (50 clusters with eps=1h)
- Outlier Detection: 0.1-4.6 seconds (IQR: 0.1-0.2s, Z-score: 2.9-4.6s)
- Technical Aggregation: 0.3-81.5 seconds (Median: 0.3-0.6s, Mean: 80-81.5s)
- **Total Pipeline**: 8-16 seconds (median aggregation)

**project_soy_2024-05 (2,064 rows)**
- Total processing: 2-3 seconds (~8 MB memory)

## Data Processing Pipeline

1. **Data Preprocessing**: Removes single-level factors, ensures consistent metadata types
2. **Logit Transformation**: Converts [0,1] range to (-∞, +∞) for percentage variables
3. **DBSCAN Clustering**: Groups timestamps, creates `dbscan_cluster` factor
4. **Outlier Detection**: Per-cell detection within factor groups
   - IQR method: Q1 - 1.5×IQR to Q3 + 1.5×IQR
   - Z-score method: mean ± 2.5×SD
5. **Technical Aggregation**: Combines multiple measurements per unit per timepoint
6. **Cluster Removal**: Removes entire DBSCAN clusters identified as outliers

## Output Files

**Processed Data Files** (saved in `data/project_NAME/`)
- `project_NAME_merged_table.rds`: Main processed dataset
- `project_NAME_vector_of_groups.rds`: Temporal grouping information
- `project_NAME_config.json`: Processing configuration

## Troubleshooting

**Validation Failures**
- Missing files: Ensure all required files present
- File format errors: Check CSV encoding and Excel format
- Column name mismatches: Verify required column names

**Processing Errors**
- DBSCAN failures: Adjust epsilon parameter or check timestamp format
- Outlier detection issues: Verify factor column names and data types
- Memory issues: Use median aggregation for large datasets

**Log Files**
- `logs/app.log`: General application logs
- `logs/error.log`: Error details
- `logs/benchmark.csv`: Performance metrics

---

**Next Steps**: After processing, proceed to the [Main Application](app_guide.md) for statistical analysis, or review [Statistical Methods](common_guide.md) for technical details.