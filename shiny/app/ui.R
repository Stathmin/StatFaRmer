# StatFaRmer User Interface
# Iteration 3: Modular Shiny

# Load required libraries
library(shiny)
library(shinyWidgets)
library(DT)
library(bslib)

# Load here package for proper paths
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here')
}
library(here)

# Load utility functions
source(here('shiny', 'common', 'utils.R'))

# Load new modules
source(here('shiny', 'app', 'module_effect_sizes.R'))
source(here('shiny', 'app', 'module_growth_summaries.R'))

# Get available projects
available_projects <- getAvailableProjects()

# Default project honoring environment override and statfarmer.project option
default_project <- if (!is.null(getOption('statfarmer.project'))) {
  getOption('statfarmer.project')
} else {
  getDefaultProject(available_projects)
}

# Load default project data for UI initialization
data_list <- tryCatch({
  loadProjectData(default_project)
}, error = function(e) {
  # Fallback to old method if project-specific files don't exist
  loadShinyData()
})

merged_table <- data_list$merged_table
vector_of_groups <- data_list$vector_of_groups

# Prepare real UI selections
ui_selections <- prepareUISelections(merged_table)

# =============================================================================
# USER INTERFACE DEFINITION
# =============================================================================

ui <- fluidPage(
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
  
  # Application title
  titlePanel("StatFaRmer"),
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      # Project selection
      selectInput(
        'selected_project',
        'Select Project:',
        choices = available_projects,
        selected = default_project
      ),
      
      # ANOVA factors selection
      selectInput(
        'anova_factors',
        'ANOVA factors:',
        choices = ui_selections$anova_factors,
        multiple = TRUE,
        selected = c('treatment', 'dbscan_cluster')
      ),
      
      # Tukey factors selection
      selectInput(
        'tukey_factors',
        'Tukey factors:',
        choices = ui_selections$anova_factors,
        multiple = TRUE,
        selected = c('treatment', 'dbscan_cluster')
      ),

      # Color by factor selection
      selectInput(
        'color_by',
        'Color by:',
        choices = c(ui_selections$anova_factors, 'Group Letter' = 'group_letter'),
        selected = 'treatment'
      ),
      
      # Grouping factor selection
      selectInput(
        'factor_grouping',
        'Grouping factor:',
        choices = ui_selections$anova_factors,
        selected = 'treatment'
      ),
      
      # Factor levels selection
      selectizeInput(
        'factor_levels',
        'Factor levels:',
        choices = ui_selections$treatments,
        multiple = TRUE,
        selected = ui_selections$treatments
      ),
      
      # Treatments selection
      selectInput(
        'treatments',
        'Selected treatments:',
        choices = ui_selections$treatments,
        multiple = TRUE,
        selected = ui_selections$treatments
      ),
      
      # Cultivars selection (default: all)
      selectInput(
        'cultivars',
        'Selected cultivars:',
        choices = ui_selections$cultivars,
        multiple = TRUE,
        selected = ui_selections$cultivars
      ),
      
      # Time clusters selection
      shinyWidgets::pickerInput(
        inputId = "timestamp_groups",
        label = "Selected time clusters:",
        choices = ui_selections$named_timestamp_groups,
        selected = {
          n <- length(ui_selections$timestamp_groups)
          if (n >= 3) {
            # Select first, middle, last as character to match picker format
            as.character(ui_selections$timestamp_groups[c(1, ceiling(n/2), n)])
          } else {
            as.character(ui_selections$timestamp_groups)
          }
        },
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10),
        multiple = TRUE
      ),
      
      # Output variables selection
      selectInput(
        'out_variables',
        'Selected trait:',
        choices = ui_selections$out_variables,
        multiple = FALSE,
        selected = ui_selections$out_variables[1]
      ),
      
      # Facet formula input
      textInput(
        'facet_formula',
        'Facet formula:',
        value = 'treatment ~ dbscan_cluster',
        placeholder = 'treatment ~ dbscan_cluster'
      ),
      
      # Plot options
      checkboxInput('timeseries_plot', ': plot timeseries with medians', FALSE),
      checkboxInput('outliers_plot', ': plot with outliers', FALSE),
      
      # Advanced Model Options
      br(),
      tags$details(
        tags$summary(
          style = "cursor: pointer; font-weight: bold; color: #337ab7;",
          "⚙️ Advanced Model Options"
        ),
        br(),
        radioButtons(
          'model_selection_method',
          'Model Selection:',
          choices = c(
            'Auto (Recommended)' = 'auto',
            'Force ANOVA' = 'aov',
            'Force LMM' = 'lmer',
            'Force Spline LMM' = 'spline'
          ),
          selected = 'auto'
        ),
        helpText(
          tags$small(
            tags$strong("Auto:"), " ≤10 timepoints → ANOVA, >10 → Spline LMM", tags$br(),
            tags$strong("ANOVA:"), " Classic fixed-effects model", tags$br(),
            tags$strong("LMM:"), " Linear mixed model with random effects", tags$br(),
            tags$strong("Spline:"), " Nonlinear growth curves (slower)"
          )
        )
      ),
      br(),
      
      # Submit button
      actionButton("submit", "Submit"),
      
      # Plot dimensions
      numericInput('plot_width', "Plot width, mm", 180),
      numericInput('plot_height', "Plot height, mm", 112),
      
      # Download button
      downloadButton("savePlot", "Save Plot as SVG")
    ),
    
    # Main panel
    mainPanel(
      # Formula display
      uiOutput("formula"),
      
      # Main plot
      plotOutput("distPlot", height = "800px", width = '1200px'),
      
      # Tabbed results
      tabsetPanel(selected = 'Config',
        tabPanel(
          "Config",
          verbatimTextOutput("project_config_json")
        ),
        tabPanel(
          "Raw Table",
          uiOutput("RAW"),
          verbatimTextOutput("raw"),
          DT::DTOutput("raw_flex"),
          downloadButton("raw_flex_downloadData", "Download Full Results")
        ),
        tabPanel(
          "Descriptive",
          uiOutput("DESCvarNames"),
          DT::DTOutput("DescriptiveTable"),
          downloadButton("descriptiveTable_downloadData", "Download Full Results")
        ),
        tabPanel(
          "ANOVA",
          uiOutput("ANOVAvarNames"),
          uiOutput("model_info_box"),
          DT::DTOutput("anovaTable"),
          downloadButton("anovaTable_downloadData", "Download Full Results"),
          verbatimTextOutput("anovaResults"),
          plotOutput("ANOVAPlot")
        ),
        tabPanel(
          "Tukey",
          uiOutput("TUKEYvarNames"),
          DT::DTOutput("tukeyTable"),
          downloadButton("tukeyTable_downloadData", "Download Full Results")
        ),
        tabPanel(
          "Assumptions",
          DT::DTOutput("assumptionsTable")
        ),
        tabPanel(
          "Group Letters",
          uiOutput("TUKEYLvarNames"),
          DT::DTOutput("tukeyLetters"),
          downloadButton("tukeyLetters_downloadData", "Download Full Results")
        ),
        
        # New scientific modules
        effectSizesUI("effect_sizes"),
        growthSummariesUI("growth_summaries")
      )
    )
  )
)
