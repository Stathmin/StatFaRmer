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
# Shared ggplot helpers (main distPlot + BLUP/DR sidecar Plots tabs)
source(here('shiny', 'common', 'plotting.R'))

# Load new modules
source(here('shiny', 'app', 'module_effect_sizes.R'))
# Curve-bridge helpers (Growth Summaries → BLUP/DR); txy_parse before plant_collapse
source(here('shiny', 'common', 'txy_parse.R'))
source(here('shiny', 'common', 'plant_collapse.R'))
source(here('shiny', 'app', 'module_growth_summaries.R'))
source(here('shiny', 'common', 'spatial_residualize.R'))
source(here('shiny', 'common', 'blup_ranking.R'))
source(here('shiny', 'app', 'module_blup_ranking.R'))
source(here('shiny', 'common', 'trait_reduction.R'))
source(here('shiny', 'app', 'module_trait_reduction.R'))


# Get available projects
available_projects <- getAvailableProjects()

# Default project honoring environment override and statfarmer.project option
default_project <- if (!is.null(getOption('statfarmer.project'))) {
  getOption('statfarmer.project')
} else {
  getDefaultProject(available_projects)
}

# Soft-fail: missing/unprocessed project must not hard-kill Shiny startup
load_result <- tryCatch({
  if (is.na(default_project) || !nzchar(as.character(default_project))) {
    stop('No projects found under data/. Run Master Wizard to create one.')
  }
  list(data = loadProjectData(default_project), error = NULL)
}, error = function(e) {
  msg <- conditionMessage(e)
  message('Failed to load selected project data (soft-fail): ', msg)
  list(data = emptyProjectData(), error = msg)
})
data_list <- load_result$data
project_load_error <- load_result$error

merged_table <- data_list$merged_table
vector_of_groups <- data_list$vector_of_groups
project_data_ready <- isTRUE(nrow(merged_table) > 0)
project_choices <- projectSelectChoices(available_projects)

# Prepare real UI selections (empty-safe)
ui_selections <- prepareUISelections(merged_table)

# Defaults that do not assume a multi-level `treatment` column (soy-325 has none)
n_treat_ui <- length(ui_selections$treatments)
has_multi_treatment <- n_treat_ui > 1L && "treatment" %in% ui_selections$anova_factors
default_anova <- intersect(
  c(if (has_multi_treatment) "treatment", "dbscan_cluster"),
  ui_selections$anova_factors
)
if (length(default_anova) == 0L) {
  default_anova <- head(setdiff(ui_selections$anova_factors, "unit"), 2L)
}
default_color <- if (has_multi_treatment) {
  "treatment"
} else if (length(setdiff(ui_selections$anova_factors, "unit")) > 0L) {
  setdiff(ui_selections$anova_factors, "unit")[1]
} else {
  "none"
}
default_grouping <- if (n_treat_ui <= 1L && "cultivar" %in% ui_selections$anova_factors) {
  "cultivar"
} else if (has_multi_treatment) {
  "treatment"
} else if (length(setdiff(ui_selections$anova_factors, "unit")) > 0L) {
  setdiff(ui_selections$anova_factors, "unit")[1]
} else {
  "none"
}
default_facet <- if (n_treat_ui <= 1L) {
  if ("dbscan_cluster" %in% ui_selections$anova_factors) "~ dbscan_cluster" else "~ ."
} else {
  "treatment ~ dbscan_cluster"
}

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

  # Empty / unprocessed project guidance (startup soft-fail)
  if (!isTRUE(project_data_ready)) {
    div(
      class = 'alert alert-warning',
      role = 'alert',
      style = 'margin: 12px 16px;',
      tags$strong('No processed project data loaded. '),
      if (!is.null(project_load_error)) {
        tags$span(project_load_error, tags$br())
      },
      tags$span(
        'Pick a processed project in the sidebar, or run Master Wizard ',
        '(launch menu option [1]) to create ',
        tags$code('*_merged_table.rds'), ' and ', tags$code('*_vector_of_groups.rds'), '.'
      )
    )
  },
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      # ----- Data (filters) -----
      h5("Data", style = "margin-top:0;"),
      selectInput(
        'selected_project',
        'Project:',
        choices = project_choices,
        selected = if (!is.na(default_project)) default_project else NULL
      ),
      selectInput(
        'out_variables',
        'Trait:',
        choices = ui_selections$out_variables,
        multiple = FALSE,
        selected = ui_selections$out_variables[1]
      ),
      shinyWidgets::pickerInput(
        inputId = "timestamp_groups",
        label = "Time clusters:",
        choices = ui_selections$named_timestamp_groups,
        selected = {
          n <- length(ui_selections$timestamp_groups)
          if (n >= 3) {
            as.character(ui_selections$timestamp_groups[c(1, ceiling(n/2), n)])
          } else {
            as.character(ui_selections$timestamp_groups)
          }
        },
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10),
        multiple = TRUE
      ),
      selectInput(
        'treatments',
        'Treatments (filter):',
        choices = ui_selections$treatments,
        multiple = TRUE,
        selected = ui_selections$treatments
      ),
      selectInput(
        'cultivars',
        'Cultivars (filter):',
        choices = ui_selections$cultivars,
        multiple = TRUE,
        selected = ui_selections$cultivars
      ),
      checkboxInput('outliers_plot', 'Include outliers in data/plot', FALSE),

      hr(),
      # ----- Model (Submit runs analysis) -----
      h5("Model"),
      helpText(tags$small("These drive ANOVA/Tukey on Submit — not plot layout.")),
      selectInput(
        'anova_factors',
        'ANOVA factors:',
        choices = ui_selections$anova_factors,
        multiple = TRUE,
        selected = default_anova
      ),
      selectInput(
        'tukey_factors',
        'Tukey / letters factors:',
        choices = ui_selections$anova_factors,
        multiple = TRUE,
        selected = default_anova
      ),
      tags$details(
        tags$summary(
          style = "cursor: pointer; font-weight: bold; color: #337ab7;",
          "Advanced model options"
        ),
        br(),
        radioButtons(
          'model_selection_method',
          'Model selection:',
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
            tags$strong("Unfit designs"), " (no multi-level factors) are skipped — no false ANOVA."
          )
        )
      ),
      actionButton("submit", "Run analysis", class = "btn-primary", width = "100%"),
      uiOutput("analysis_status_banner"),

      hr(),
      # ----- Plot (visuals only) -----
      h5("Plot"),
      helpText(tags$small("Visuals only. CLD letters appear after a successful Tukey run.")),
      selectInput(
        'color_by',
        'Color by:',
        choices = c(ui_selections$anova_factors, 'Group Letter' = 'group_letter', 'None' = 'none'),
        selected = default_color
      ),
      selectInput(
        'factor_grouping',
        'X-axis grouping:',
        choices = c('None' = 'none', ui_selections$anova_factors),
        selected = default_grouping
      ),
      selectizeInput(
        'factor_levels',
        'X-axis levels (optional):',
        choices = ui_selections$treatments,
        multiple = TRUE,
        selected = ui_selections$treatments
      ),
      textInput(
        'facet_formula',
        'Facet formula:',
        value = default_facet,
        placeholder = '~ dbscan_cluster  |  treatment ~ cultivar'
      ),
      numericInput(
        'facet_ncol',
        'Facet wrap columns:',
        value = 2, min = 1, max = 8, step = 1
      ),
      checkboxInput('timeseries_plot', 'Timeseries (medians)', FALSE),
      numericInput('plot_width', "Export width, mm", 180),
      numericInput('plot_height', "Export height, mm", 112),
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
        growthSummariesUI("growth_summaries"),
        blupRankingUI("blup_ranking"),
        traitReductionUI("trait_reduction")
      )
    )
  )
)
