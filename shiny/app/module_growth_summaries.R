# Growth Summaries Module
# Calculates AUC, peak, growth_rate, last (and optional mean) over selected times.
# Multi-trait panel → trait__metric columns for Trait DR / BLUP (soy-325 style).

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
        p("Collapse selected times to named plant-level metrics (trait__auc / __peak / __growth_rate / __last) for Trait DR and BLUP."),
        tags$div(
          class = "alert alert-info",
          role = "alert",
          strong("Path to multi-trait PCA / BLUP: "),
          "Pick several PlantEye traits below → Calculate Summaries → Trait DR loadings span those traits (not just one sidebar curve)."
        ),
        hr()
      )
    ),

    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          p(strong("Focus trait (plots):"), textOutput(ns("current_trait"), inline = TRUE)),
          br(),
          selectInput(
            ns("panel_traits"),
            "Traits for plant panel (Trait DR / BLUP):",
            choices = c(),
            multiple = TRUE,
            selectize = TRUE
          ),
          helpText(tags$small(
            "Each selected trait becomes columns like digital_biomass_mm3__auc. ",
            "Focus trait (left bar) is always included. Same basis as soy-325 multi-trait PCA."
          )),
          selectInput(ns("growth_grouping"), "Group by:",
                     choices = c(), multiple = TRUE),
          hr(),
          h5("Metrics to calculate:"),
          checkboxInput(ns("calc_auc"), "Area Under Curve (AUC)", value = TRUE),
          checkboxInput(ns("calc_peak"), "Peak value", value = TRUE),
          checkboxInput(ns("calc_rate"), "Growth rate (Δy / Δdays)", value = TRUE),
          checkboxInput(ns("calc_last"), "Endpoint (last selected time)", value = TRUE),
          checkboxInput(ns("calc_mean"), "Mean over selected times", value = FALSE),
          helpText(
            tags$small(
              tags$strong("Selected times only:"), " collapse uses the main timestamp filter.", tags$br(),
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
            p("One row per plant; columns are trait__metric. Download CSV or continue to Trait DR / BLUP."),
            uiOutput(ns("summary_note")),
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
              tags$li(strong("AUC:"), " Trapezoidal area under each plant’s trait-vs-time curve (selected times; time in hours)."),
              tags$li(strong("Peak:"), " Maximum value in the selected window."),
              tags$li(strong("Growth rate:"), " (last − first) / Δdays over selected times."),
              tags$li(strong("Last:"), " Value at the last selected timestamp."),
              tags$li(strong("Mean:"), " Mean of selected time points (explicit metric only — not silent DR default).")
            ),
            br(),
            h5("Why multi-trait?"),
            p("Trait DR PCA needs several phenotype columns. Collapsing many traits to trait__metric ",
              "matches the offline soy-325 suite so loadings can mix biomass, NDVI, hull, GLI, etc.")
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

    output$current_trait <- renderText({
      req(combined_inputs$out_variables)
      combined_inputs$out_variables
    })

    observe({
      req(combined_inputs$anova_factors)
      choices <- setdiff(combined_inputs$anova_factors, "dbscan_cluster")
      preferred <- intersect(c("treatment", "cultivar", "genotype"), choices)
      selected <- if (length(preferred) > 0L) preferred else head(choices, 1L)
      updateSelectInput(session, "growth_grouping",
                       choices = choices,
                       selected = selected)
    })

    # Populate multi-trait panel from current filtered / project table
    observe({
      data <- tryCatch(filteredData(), error = function(e) NULL)
      if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
        pd <- tryCatch(projectData(), error = function(e) NULL)
        if (!is.null(pd) && is.data.frame(pd$merged_table)) {
          data <- pd$merged_table
        }
      }
      cands <- .sf_gs_phenotype_candidates(data)
      focus <- tryCatch(combined_inputs$out_variables, error = function(e) NULL)
      if (length(cands) < 1L) {
        updateSelectInput(session, "panel_traits", choices = character(0), selected = character(0))
        return()
      }
      default <- .sf_gs_default_panel_traits(cands, focus_trait = focus, max_n = 12L)
      cur <- isolate(input$panel_traits)
      keep <- intersect(cur, cands)
      selected <- if (length(keep) >= 1L) {
        unique(c(intersect(as.character(focus), cands), keep))
      } else {
        default
      }
      updateSelectInput(session, "panel_traits", choices = cands, selected = selected)
    })

    growth_summaries <- eventReactive(input$calculate_growth, {
      req(filteredData())
      req(nrow(filteredData()) > 0)

      focus <- combined_inputs$out_variables
      data <- filteredData()

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

      traits <- unique(as.character(input$panel_traits))
      traits <- intersect(traits, names(data))
      if (!is.null(focus) && nzchar(as.character(focus)[[1L]]) &&
          as.character(focus)[[1L]] %in% names(data)) {
        traits <- unique(c(as.character(focus)[[1L]], traits))
      }
      if (length(traits) < 1L) {
        traits <- intersect(.sf_gs_phenotype_candidates(data), names(data))
        traits <- head(traits, 8L)
      }
      shiny::validate(shiny::need(
        length(traits) >= 1L,
        "Select at least one phenotype for the plant panel."
      ))

      metrics <- character(0)
      if (isTRUE(input$calc_auc)) metrics <- c(metrics, "auc")
      if (isTRUE(input$calc_peak)) metrics <- c(metrics, "peak")
      if (isTRUE(input$calc_rate)) metrics <- c(metrics, "growth_rate")
      if (isTRUE(input$calc_last)) metrics <- c(metrics, "last")
      if (isTRUE(input$calc_mean)) metrics <- c(metrics, "mean")
      shiny::validate(shiny::need(
        length(metrics) >= 1L,
        "Enable at least one metric (AUC / peak / growth_rate / last / mean)."
      ))

      withProgress(message = "Collapsing multi-trait plant panel…", value = 0.15, {
        incProgress(
          0.25,
          detail = sprintf("%d traits × %s | %d time cluster(s)",
                           length(traits), paste(metrics, collapse = "+"),
                           length(selected_clusters))
        )

        summaries <- collapseManyTraitsBySelectedTimes(
          data = data,
          traits = traits,
          time_col = "dbscan_cluster",
          time_values = selected_clusters,
          metrics = metrics
        )

        attr(summaries, "focus_trait") <- if (!is.null(focus)) as.character(focus)[[1L]] else NA_character_
        attr(summaries, "panel_traits") <- traits
        attr(summaries, "gs_metrics") <- metrics
        attr(summaries, "selected_clusters") <- selected_clusters

        # Optional n_obs for focus trait series
        if ("unit" %in% names(data) && "unit" %in% names(summaries)) {
          series <- data[data$dbscan_cluster %in% selected_clusters, , drop = FALSE]
          n_tbl <- as.data.frame(table(series$unit), stringsAsFactors = FALSE)
          names(n_tbl) <- c("unit", "n_obs")
          n_tbl$unit <- as.character(n_tbl$unit)
          summaries$unit <- as.character(summaries$unit)
          summaries <- merge(summaries, n_tbl, by = "unit", all.x = TRUE)
        }

        incProgress(0.5, detail = paste(nrow(summaries), "plants |",
                                        length(attr(summaries, "metric_cols")), "metric cols"))
        summaries
      })
    })

    output$summary_note <- renderUI({
      if (!isTruthy(input$calculate_growth) || input$calculate_growth < 1L) {
        return(NULL)
      }
      gs <- tryCatch(growth_summaries(), error = function(e) NULL)
      if (is.null(gs)) return(NULL)
      mc <- attr(gs, "metric_cols")
      tr <- attr(gs, "panel_traits")
      tags$div(
        class = "alert alert-success",
        role = "alert",
        sprintf(
          "%d plants · %d traits · %d named metrics (e.g. %s). Trait DR can run multi-trait PCA on these columns.",
          nrow(gs),
          length(tr %||% character(0)),
          length(mc %||% character(0)),
          paste(utils::head(mc %||% character(0), 3L), collapse = ", ")
        )
      )
    })

    output$growth_summary_table <- DT::renderDataTable({
      if (!isTruthy(input$calculate_growth) || input$calculate_growth < 1L) {
        return(DT::datatable(
          data.frame(
            Message = "Click Calculate Summaries to build trait__metric plant rows (required before Trait DR / BLUP)."
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

    .sf_gs_plot_df <- function(metric) {
      data <- growth_summaries()
      focus <- attr(data, "focus_trait")
      if (is.null(focus) || !nzchar(as.character(focus)[[1L]])) {
        focus <- combined_inputs$out_variables
      }
      ycol <- .sf_gs_metric_column(data, metric, focus_trait = focus)
      list(data = data, ycol = ycol, focus = focus)
    }

    output$auc_plot <- renderPlot({
      req(growth_summaries(), input$calc_auc)
      pack <- .sf_gs_plot_df("auc")
      if (is.na(pack$ycol)) {
        plot.new()
        text(0.5, 0.5, "No AUC column — enable AUC and recalculate", cex = 1.3)
        return()
      }
      data <- pack$data
      grp <- input$growth_grouping
      grp <- grp[!is.null(grp) & !is.na(grp) & nzchar(as.character(grp))]
      grp <- intersect(as.character(grp), names(data))
      if (length(grp) < 1L) {
        plot.new()
        text(0.5, 0.5, "Select at least one grouping factor present in the data", cex = 1.2)
        return()
      }
      ycol <- pack$ycol
      p <- ggplot(data, aes(x = !!sym(grp[1]), y = !!sym(ycol))) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = if (length(grp) >= 2) NA else "#4CAF50") +
        geom_jitter(width = 0.2, alpha = 0.6)
      if (length(grp) >= 2) {
        p <- ggplot(data, aes(x = !!sym(grp[1]), y = !!sym(ycol), color = !!sym(grp[2]))) +
          geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6)
      }
      if (length(grp) >= 3) {
        p <- p + facet_wrap(as.formula(paste("~", grp[3])))
      }
      p +
        getStatfarmerTheme() +
        labs(
          title = paste0("AUC by group — ", ycol),
          subtitle = "Selected times only (trapezoidal, hours)",
          x = grp[1],
          y = ycol
        ) +
        theme(text = element_text(size = 14))
    })

    output$rate_plot <- renderPlot({
      req(growth_summaries(), input$calc_rate)
      pack <- .sf_gs_plot_df("growth_rate")
      if (is.na(pack$ycol)) {
        plot.new()
        text(0.5, 0.5, "No growth_rate column — enable growth rate and recalculate", cex = 1.3)
        return()
      }
      data <- pack$data
      grp <- input$growth_grouping
      grp <- grp[!is.null(grp) & !is.na(grp) & nzchar(as.character(grp))]
      grp <- intersect(as.character(grp), names(data))
      if (length(grp) < 1L) {
        plot.new()
        text(0.5, 0.5, "Select at least one grouping factor present in the data", cex = 1.2)
        return()
      }
      ycol <- pack$ycol
      p <- ggplot(data, aes(x = !!sym(grp[1]), y = !!sym(ycol))) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, fill = if (length(grp) >= 2) NA else "#2196F3") +
        geom_jitter(width = 0.2, alpha = 0.6)
      if (length(grp) >= 2) {
        p <- ggplot(data, aes(x = !!sym(grp[1]), y = !!sym(ycol), color = !!sym(grp[2]))) +
          geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.6)
      }
      if (length(grp) >= 3) {
        p <- p + facet_wrap(as.formula(paste("~", grp[3])))
      }
      p +
        getStatfarmerTheme() +
        labs(
          title = paste0("Growth rate by group — ", ycol),
          subtitle = "Selected times only (Δy / Δdays)",
          x = grp[1],
          y = ycol
        ) +
        theme(text = element_text(size = 14))
    })

    output$download_summaries <- downloadHandler(
      filename = function() {
        paste0("growth_summaries_multitrait_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(growth_summaries(), file, row.names = FALSE)
      }
    )

    list(
      plant_summaries = growth_summaries,
      growth_summaries = growth_summaries
    )
  })
}
