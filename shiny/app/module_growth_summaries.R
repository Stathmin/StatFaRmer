# Growth Summaries Module
# Calculates AUC, time-to-peak, growth rates, and other temporal metrics

#' Growth Summaries UI
#' @param id Module ID
growthSummariesUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Growth Summaries",
    fluidRow(
      column(12,
        h3("Growth Trajectory Summaries"),
        p("Summarize growth dynamics: AUC, peak values, growth rates, and timing."),
        hr()
      )
    ),
    
    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          p(strong("Growth variable:"), textOutput(ns("current_trait"), inline = TRUE)),
          br(),
          selectInput(ns("growth_grouping"), "Group by:", 
                     choices = c(), multiple = TRUE),
          hr(),
          h5("Metrics to calculate:"),
          checkboxInput(ns("calc_auc"), "Area Under Curve (AUC)", value = TRUE),
          checkboxInput(ns("calc_peak"), "Peak value & timing", value = TRUE),
          checkboxInput(ns("calc_rate"), "Growth rate (slope)", value = TRUE),
          checkboxInput(ns("calc_relative"), "Relative growth rate", value = FALSE),
          hr(),
          actionButton(ns("calculate_growth"), "Calculate Summaries", 
                      class = "btn-primary", icon = icon("chart-line"))
        )
      ),
      
      column(8,
        tabsetPanel(
          tabPanel("Summary Table",
            br(),
            DT::dataTableOutput(ns("growth_summary_table")),
            br(),
            downloadButton(ns("download_summaries"), "Download CSV")
          ),
          tabPanel("AUC Comparison",
            br(),
            plotOutput(ns("auc_plot"), height = "500px")
          ),
          tabPanel("Growth Rates",
            br(),
            plotOutput(ns("rate_plot"), height = "500px")
          ),
          tabPanel("About",
            br(),
            h4("Growth Summary Metrics"),
            tags$ul(
              tags$li(strong("AUC (Area Under Curve):"), " Total accumulated growth over time. Higher AUC = more total growth."),
              tags$li(strong("Peak value:"), " Maximum value reached during observation period."),
              tags$li(strong("Time to peak:"), " When maximum value was reached (early/late grower)."),
              tags$li(strong("Growth rate (slope):"), " Linear rate of increase (units per time)."),
              tags$li(strong("Relative growth rate (RGR):"), " Proportional growth rate: (ln(final) - ln(initial)) / time.")
            ),
            br(),
            h5("Why use growth summaries?"),
            p("Growth summaries reduce time-series data to interpretable metrics. This allows:"),
            tags$ul(
              tags$li("Fast ANOVA on summarized values (no time factor)"),
              tags$li("Direct interpretation: 'Treatment X increased total growth by 25%'"),
              tags$li("Avoid complex mixed models for simple questions")
            )
          )
        )
      )
    )
  )
}

#' Growth Summaries Server
#' @param id Module ID
#' @param filteredData Reactive filtered data
#' @param combined_inputs Reactive combined inputs
#' @param projectData Reactive project data  
growthSummariesServer <- function(id, filteredData, combined_inputs, projectData) {
  moduleServer(id, function(input, output, session) {
    
    # Display current trait
    output$current_trait <- renderText({
      req(combined_inputs$out_variables)
      combined_inputs$out_variables
    })
    
    # Update choices
    observe({
      req(combined_inputs$anova_factors)
      updateSelectInput(session, "growth_grouping",
                       choices = setdiff(combined_inputs$anova_factors, "dbscan_cluster"),
                       selected = c("treatment", "cultivar"))
    })
    
    # Calculate growth summaries
    growth_summaries <- eventReactive(input$calculate_growth, {
      req(projectData(), combined_inputs$out_variables)
      
      # Use full data (not filtered) to get complete time series per unit
      data <- projectData()$merged_table
      metric <- combined_inputs$out_variables  # Use selected trait from main UI
      grouping <- input$growth_grouping
      
      # Ensure dbscan_cluster is sorted by time
      data <- data %>% arrange(timestamp_group)
      
      # Calculate per unit
      summaries <- data %>%
        group_by(unit, across(all_of(grouping))) %>%
        summarise(
          n_obs = n(),
          first_value = first(!!sym(metric)),
          final_value = last(!!sym(metric)),
          .groups = "drop"
        )
      
      # AUC (trapezoidal rule)
      if (input$calc_auc) {
        auc_data <- data %>%
          arrange(unit, timestamp_group) %>%
          group_by(unit) %>%
          summarise(
            auc = sum(diff(as.numeric(timestamp_group)) / 3600 * 
                     (head(!!sym(metric), -1) + tail(!!sym(metric), -1)) / 2,
                     na.rm = TRUE),
            .groups = "drop"
          )
        summaries <- summaries %>% left_join(auc_data, by = "unit")
      }
      
      # Peak value and timing
      if (input$calc_peak) {
        peak_data <- data %>%
          group_by(unit) %>%
          slice_max(!!sym(metric), n = 1, with_ties = FALSE) %>%
          summarise(
            peak_value = !!sym(metric),
            time_to_peak = as.numeric(difftime(timestamp_group, min(data$timestamp_group), units = "days")),
            .groups = "drop"
          )
        summaries <- summaries %>% left_join(peak_data, by = "unit")
      }
      
      # Growth rate (linear slope)
      if (input$calc_rate) {
        rate_data <- data %>%
          group_by(unit) %>%
          summarise(
            time_range = as.numeric(difftime(max(timestamp_group), min(timestamp_group), units = "days")),
            value_change = last(!!sym(metric)) - first(!!sym(metric)),
            growth_rate = value_change / time_range,
            .groups = "drop"
          )
        summaries <- summaries %>% left_join(rate_data[, c("unit", "growth_rate")], by = "unit")
      }
      
      # Relative growth rate
      if (input$calc_relative) {
        summaries <- summaries %>%
          mutate(
            rgr = (log(final_value) - log(first_value)) / n_obs
          )
      }
      
      summaries
    })
    
    # Render summary table
    output$growth_summary_table <- DT::renderDataTable({
      req(growth_summaries())
      DT::datatable(
        growth_summaries() %>%
          mutate(across(where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    # AUC plot
    output$auc_plot <- renderPlot({
      req(growth_summaries(), input$calc_auc)
      
      if (!"auc" %in% names(growth_summaries())) {
        plot.new()
        text(0.5, 0.5, "Enable 'Area Under Curve' to see this plot", cex = 1.5)
        return()
      }

      grp <- input$growth_grouping
      data <- growth_summaries()

      # Base plot on first grouping
      p <- ggplot(data, aes(x = !!sym(grp[1]), y = auc)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = if (length(grp) >= 2) NA else "#4CAF50") +
        geom_jitter(width = 0.2, alpha = 0.6)

      # If second grouping provided, color by it
      if (length(grp) >= 2) {
        p <- ggplot(data, aes(x = !!sym(grp[1]), y = auc, color = !!sym(grp[2]))) +
          geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6)
      }

      # If third grouping provided, facet by it
      if (length(grp) >= 3) {
        p <- p + facet_wrap(as.formula(paste("~", grp[3])))
      }

      p +
        getStatfarmerTheme() +
        labs(title = "Area Under Curve by Group",
             x = grp[1],
             y = "AUC") +
        theme(text = element_text(size = 14))
    })
    
    # Growth rate plot
    output$rate_plot <- renderPlot({
      req(growth_summaries(), input$calc_rate)
      
      if (!"growth_rate" %in% names(growth_summaries())) {
        plot.new()
        text(0.5, 0.5, "Enable 'Growth rate (slope)' to see this plot", cex = 1.5)
        return()
      }

      grp <- input$growth_grouping
      data <- growth_summaries()

      p <- ggplot(data, aes(x = !!sym(grp[1]), y = growth_rate)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = if (length(grp) >= 2) NA else "#2196F3") +
        geom_jitter(width = 0.2, alpha = 0.6)

      if (length(grp) >= 2) {
        p <- ggplot(data, aes(x = !!sym(grp[1]), y = growth_rate, color = !!sym(grp[2]))) +
          geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6)
      }

      if (length(grp) >= 3) {
        p <- p + facet_wrap(as.formula(paste("~", grp[3])))
      }

      p +
        getStatfarmerTheme() +
        labs(title = "Growth Rate by Group",
             x = grp[1],
             y = "Growth Rate (units/day)") +
        theme(text = element_text(size = 14))
    })
    
    # Download handler
    output$download_summaries <- downloadHandler(
      filename = function() {
        paste0("growth_summaries_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(growth_summaries(), file, row.names = FALSE)
      }
    )
    
  })
}

