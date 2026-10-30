# Growth Summaries Module
# Calculates AUC, time-to-peak, growth rates, and other temporal metrics
# Curve bridge: metrics use user-selected time clusters only (not full series).

#' Map combined_inputs timestamp selection to dbscan_cluster values.
#' @param data Plant × time data.frame (e.g. merged_table)
#' @param timestamp_groups Selected timestamp_group values (POSIXct or character)
#' @return Unique dbscan_cluster values for the selection (may be empty)
.sf_selected_dbscan_clusters <- function(data, timestamp_groups) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
    return(unique(data$dbscan_cluster)[0])
  }
  if (!all(c("timestamp_group", "dbscan_cluster") %in% names(data))) {
    if ("dbscan_cluster" %in% names(data)) {
      return(unique(data$dbscan_cluster))
    }
    return(NULL)
  }
  if (is.null(timestamp_groups) || length(timestamp_groups) < 1L) {
    return(unique(data$dbscan_cluster)[0])
  }
  tg <- timestamp_groups
  if (is.character(tg) && inherits(data$timestamp_group, "POSIXct")) {
    tg <- as.POSIXct(tg, tz = "UTC")
  }
  unique(data$dbscan_cluster[data$timestamp_group %in% tg])
}

#' Growth Summaries UI
#' @param id Module ID
growthSummariesUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Growth Summaries",
    fluidRow(
      column(12,
        h3("Growth Trajectory Summaries"),
        p("Summarize growth dynamics over the currently selected time points: AUC, peak, growth rates, and timing."),
        tags$div(
          class = "alert alert-info",
          role = "alert",
          strong("Path to BLUP / Trait DR: "),
          "Click Calculate Summaries to build one row per plant, then open Trait DR (optional PCA) or BLUP to rank genotypes. Raw multi-time rows are blocked on those tabs."
        ),
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
          helpText(
            tags$small(
              tags$strong("Selected times only:"), " collapse uses the main timestamp filter, not the full series.", tags$br(),
              tags$strong("Next steps:"), " Trait DR → PC hyper-traits, then BLUP → genotype ranks."
            )
          ),
          hr(),
          actionButton(ns("calculate_growth"), "Calculate Summaries", 
                      class = "btn-primary", icon = icon("chart-line"))
        )
      ),
      
      column(8,
        tabsetPanel(
          tabPanel("Summary Table",
            br(),
            p("One row per plant over the selected time clusters. Download CSV or continue to Trait DR / BLUP."),
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
              tags$li(strong("AUC (Area Under Curve):"), " Total accumulated growth over selected times. Higher AUC = more total growth."),
              tags$li(strong("Peak value:"), " Maximum value reached during the selected observation period."),
              tags$li(strong("Time to peak:"), " When maximum value was reached (early/late grower)."),
              tags$li(strong("Growth rate (slope):"), " Linear rate of increase (units per time)."),
              tags$li(strong("Relative growth rate (RGR):"), " Proportional growth rate: (ln(final) - ln(initial)) / time.")
            ),
            br(),
            h5("Why use growth summaries?"),
            p("Growth summaries reduce time-series data to interpretable metrics over the selected time clusters. This allows:"),
            tags$ul(
              tags$li("Fast ANOVA on summarized values (no time factor)"),
              tags$li("Direct interpretation: 'Treatment X increased total growth by 25%'"),
              tags$li("Plant-level curve traits for BLUP / Trait DR (one row per plant)"),
              tags$li("Required first step before Calculate BLUP or Calculate PCA")
            )
          )
        )
      )
    )
  )
}

#' Growth Summaries Server
#' @param id Module ID
#' @param filteredData Reactive filtered data (includes selected times)
#' @param combined_inputs Reactive combined inputs
#' @param projectData Reactive project data
#' @return List with `plant_summaries` / `growth_summaries` reactives (one row per plant)
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
    
    # Calculate growth summaries over selected times only (curve bridge)
    growth_summaries <- eventReactive(input$calculate_growth, {
      req(filteredData(), combined_inputs$out_variables)
      req(nrow(filteredData()) > 0)

      metric <- combined_inputs$out_variables
      data <- filteredData()

      # Explicit selected clusters from main timestamp filter (same as analysis)
      selected_clusters <- NULL
      if (!is.null(projectData()) && nrow(projectData()$merged_table) > 0) {
        selected_clusters <- .sf_selected_dbscan_clusters(
          projectData()$merged_table,
          combined_inputs$timestamp_groups
        )
      }
      if (is.null(selected_clusters) || length(selected_clusters) < 1L) {
        selected_clusters <- unique(data$dbscan_cluster)
      }
      req(length(selected_clusters) >= 1L)

      withProgress(message = "Collapsing plants over selected times…", value = 0.15, {
        metrics <- "mean"
        if (isTRUE(input$calc_auc)) {
          metrics <- c(metrics, "auc")
        }
        if (isTRUE(input$calc_peak)) {
          metrics <- c(metrics, "peak")
        }

        incProgress(0.25, detail = paste(length(selected_clusters), "time cluster(s)"))

        # Collapse via shared helper — selected times only, not full merged_table
        summaries <- collapsePlantsBySelectedTimes(
          data = data,
          trait = metric,
          time_col = "dbscan_cluster",
          time_values = selected_clusters,
          metrics = metrics
        )

        # Series restricted to selected clusters for rate / timing extras
        series <- data %>%
          dplyr::filter(dbscan_cluster %in% selected_clusters) %>%
          dplyr::arrange(timestamp_group)

        base_stats <- series %>%
          dplyr::group_by(unit) %>%
          dplyr::summarise(
            n_obs = dplyr::n(),
            first_value = dplyr::first(.data[[metric]]),
            final_value = dplyr::last(.data[[metric]]),
            .groups = "drop"
          )
        summaries <- summaries %>% dplyr::left_join(base_stats, by = "unit")

        if (isTRUE(input$calc_peak) && "timestamp_group" %in% names(series)) {
          t0 <- min(series$timestamp_group, na.rm = TRUE)
          peak_timing <- series %>%
            dplyr::group_by(unit) %>%
            dplyr::slice_max(.data[[metric]], n = 1, with_ties = FALSE) %>%
            dplyr::summarise(
              time_to_peak = as.numeric(difftime(timestamp_group, t0, units = "days")),
              .groups = "drop"
            )
          summaries <- summaries %>% dplyr::left_join(peak_timing, by = "unit")
        }

        if (isTRUE(input$calc_rate)) {
          rate_data <- series %>%
            dplyr::group_by(unit) %>%
            dplyr::summarise(
              time_range = as.numeric(difftime(
                max(timestamp_group), min(timestamp_group),
                units = "days"
              )),
              value_change = dplyr::last(.data[[metric]]) - dplyr::first(.data[[metric]]),
              growth_rate = value_change / time_range,
              .groups = "drop"
            )
          summaries <- summaries %>%
            dplyr::left_join(rate_data[, c("unit", "growth_rate")], by = "unit")
        }

        if (isTRUE(input$calc_relative)) {
          summaries <- summaries %>%
            dplyr::mutate(
              rgr = (log(final_value) - log(first_value)) / n_obs
            )
        }

        incProgress(0.5, detail = paste(nrow(summaries), "plants"))
        summaries
      })
    })
    
    # Render summary table
    output$growth_summary_table <- DT::renderDataTable({
      if (!isTruthy(input$calculate_growth) || input$calculate_growth < 1L) {
        return(DT::datatable(
          data.frame(
            Message = "Click Calculate Summaries to collapse plants (required before Trait DR / BLUP)."
          ),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
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
        labs(title = "Area Under Curve by Group (selected times)",
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
        labs(title = "Growth Rate by Group (selected times)",
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

    # Expose plant-level summaries for BLUP / DR (prefer return over session$userData)
    list(
      plant_summaries = growth_summaries,
      growth_summaries = growth_summaries
    )
  })
}
