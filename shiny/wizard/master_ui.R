library(shiny)
library(shinyjs)
library(plotly)
library(bslib)

master_ui <- fluidPage(
  # Professional scientific theme
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#3498db",      # Scientific blue
    secondary = "#95a5a6",    # Neutral gray
    success = "#27ae60",      # Green for success
    info = "#3498db",         # Blue for information
    warning = "#f39c12",      # Orange for warnings
    danger = "#e74c3c",       # Red for errors
    base_font = font_google("Inter"),
    code_font = font_google("Fira Code")
  ),
  
  useShinyjs(),
  titlePanel('StatFaRmer — Master Wizard'),
  sidebarLayout(
    sidebarPanel(
      h4('1) Input folder'),
      selectizeInput('input_dir', 'Project folder under data/', choices = NULL, options = list(placeholder = 'Select a project folder')), 
      actionButton('validate_btn', 'Validate input'),
      hr(),
      conditionalPanel(
        condition = 'output.validated == "true"',
        h4('2) Processing parameters'),
        numericInput('eps_hours', 'DBSCAN eps (hours)', value = 1, min = 0.1, step = 0.1),
        checkboxInput('use_logit', 'Apply logit transform', value = FALSE),
        checkboxInput('logit_treat_inf', 'Treat ±Inf/NaN (mask to NA)', value = TRUE),
        checkboxInput('use_iqr', 'Enable outlier handling by cells', value = FALSE),
        radioButtons('outlier_method', 'Outlier handling method', 
                    choices = c('Remove outliers' = 'remove', 'Winsorize outliers' = 'winsorize', 'Keep all outliers' = 'keep'), 
                    selected = 'remove'),
        radioButtons('outlier_detection', 'Detection method', 
                    choices = c('IQR' = 'iqr', 'Z-score' = 'zscore'), 
                    selected = 'iqr'),
        selectizeInput('iqr_factors', 'Factor columns for cells', choices = NULL, multiple = TRUE),
        selectizeInput('dbscan_numeric', 'Numeric variable for DBSCAN preview (y vs timestamp)', choices = NULL, multiple = FALSE),
        textInput('facet_formula', 'Facet formula (e.g., treatment ~ dbscan_cluster)', value = '~ .'),
        radioButtons('tech_agg', 'Aggregation of technical replicates', choices = c('median','mean'), selected = 'median'),
        selectizeInput('outlier_clusters', 'Select DBSCAN clusters as outliers', choices = NULL, multiple = TRUE),
        checkboxInput('hide_outliers_plot', 'Hide outliers in plot', value = FALSE),
        hr(),
        h4('3) Export'),
        textInput('export_name', 'Export name (suffix)', placeholder = 'custom'),
        actionButton('run_btn', 'Run processing and export')
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Validation', verbatimTextOutput('validation_out')),
        tabPanel('Preview params', verbatimTextOutput('params_out')),
        tabPanel('DBSCAN preview', plotly::plotlyOutput('dbscan_plot', height = '420px')),
        tabPanel('Log', verbatimTextOutput('log_tail'))
      )
    )
  )
)


