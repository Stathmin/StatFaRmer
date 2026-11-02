# Trait DR Module
# Plant × traits PCA → hyper-traits (PC scores) for BLUP response

.sf_dr_read_plants <- function(plant_summaries) {
  tryCatch({
    ps <- plant_summaries()
    if (is.null(ps) || !is.data.frame(ps) || nrow(ps) < 1L) {
      return(NULL)
    }
    ps
  }, error = function(e) NULL)
}

#' Trait DR UI
#' @param id Module ID
traitReductionUI <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Trait DR",
    fluidRow(
      column(12,
        h3("Trait Dimensionality Reduction"),
        p("PCA on plant × named Growth Summaries metrics (trait__auc / __peak / …). PC scores attach to the same plant rows for BLUP."),
        uiOutput(ns("dr_guidance")),
        hr()
      )
    ),

    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          selectInput(ns("dr_traits"), "Traits (named GS metrics):",
                      choices = c(), multiple = TRUE),
          numericInput(ns("n_comp"), "Number of PCs (blank = all):",
                       value = NA, min = 1, max = 50, step = 1),
          checkboxInput(ns("scale_traits"), "Scale traits before PCA", value = TRUE),
          helpText(
            tags$small(
              tags$strong("Rows:"), " plants from Growth Summaries (one row per plant).", tags$br(),
              tags$strong("Columns:"), " multi-trait panel metrics (e.g. digital_biomass_mm3__auc, ndvi_average__last) — not raw scans.", tags$br(),
              tags$strong("Next:"), " after Calculate PCA, pick PC1… as the trait on the BLUP tab."
            )
          ),
          hr(),
          actionButton(ns("calculate_dr"), "Calculate PCA",
                       class = "btn-primary", icon = icon("calculator"))
        )
      ),

      column(8,
        tabsetPanel(
          tabPanel("Plots",
            br(),
            p("Shared StatFaRmer ggplot theme (createDrScreePlot / Scores / Loadings in plotting.R)."),
            selectInput(ns("dr_color_by"), "Color plant scores by:",
                        choices = c("cultivar", "treatment", "none"),
                        selected = "cultivar"),
            plotOutput(ns("scree_plot"), height = "320px"),
            downloadButton(ns("scree_svg"), "Save scree SVG"),
            hr(),
            plotOutput(ns("scores_plot"), height = "400px"),
            downloadButton(ns("scores_svg"), "Save scores SVG"),
            hr(),
            plotOutput(ns("loadings_plot"), height = "420px"),
            downloadButton(ns("loadings_svg"), "Save loadings SVG")
          ),
          tabPanel("Scree table",
            br(),
            p("Proportion of variance per PC."),
            DT::dataTableOutput(ns("variance_table"))
          ),
          tabPanel("Loadings",
            br(),
            p("How each selected trait loads on each PC."),
            DT::dataTableOutput(ns("loadings_table")),
            br(),
            downloadButton(ns("loadings_download"), "Download loadings (CSV)")
          ),
          tabPanel("Plant scores",
            br(),
            p("PC columns bound to the same plant rows (unit / cultivar preserved). Use these as BLUP responses."),
            DT::dataTableOutput(ns("scores_table")),
            br(),
            downloadButton(ns("scores_download"), "Download plant scores (CSV)")
          ),
          tabPanel("Trait clusters",
            br(),
            p("Clusters from absolute correlation distance among selected traits."),
            DT::dataTableOutput(ns("cluster_table")),
            br(),
            downloadButton(ns("cluster_download"), "Download clusters (CSV)")
          ),
          tabPanel("About",
            br(),
            h4("Trait DR sidecar"),
            tags$ul(
              tags$li(strong("Rows:"), " plants only (Growth Summaries multi-trait collapse)."),
              tags$li(strong("Columns:"), " named temporal metrics trait__auc / __peak / __growth_rate / __last (and optional __mean)."),
              tags$li(strong("Plots:"), " createDrScreePlot, createDrScoresPlot, createDrLoadingsPlot via plotting.R."),
              tags$li(strong("Output:"), " PC1… on the same plant table for the BLUP trait chooser."),
              tags$li(strong("Refresh:"), " re-running Growth Summaries clears cached PCs — Calculate PCA again."),
              tags$li(strong("Anti-goal:"), " never silent across-time raw means; never genotype×trait mean matrix.")
            )
          )
        )
      )
    )
  )
}

#' Trait DR Server
#' @param id Module ID
#' @param plant_summaries Reactive plant-level summaries from growthSummariesServer
#' @param combined_inputs Reactive combined inputs (unused reserved)
#' @return List with `plant_with_pcs` reactive (NULL until Calculate PCA)
traitReductionServer <- function(id, plant_summaries, combined_inputs) {
  moduleServer(id, function(input, output, session) {

    # Respect left-sidebar cultivar/treatment filters (same bus as main plot)
    plants_sidebar <- reactive({
      ps <- .sf_dr_read_plants(plant_summaries)
      .sf_filter_plants_by_sidebar(ps, combined_inputs)
    })

    output$dr_guidance <- renderUI({
      ps_raw <- .sf_dr_read_plants(plant_summaries)
      if (is.null(ps_raw)) {
        tags$div(
          class = "alert alert-warning",
          role = "alert",
          strong("Growth Summaries required. "),
          "Open Growth Summaries and click Calculate Summaries so Trait DR has one row per plant."
        )
      } else {
        ps <- plants_sidebar()
        n_after <- if (is.null(ps)) 0L else nrow(ps)
        eligible <- .sf_dr_eligible_traits(if (!is.null(ps) && nrow(ps) > 0L) ps else ps_raw)
        if (length(eligible) < 2L) {
          tags$div(
            class = "alert alert-warning",
            role = "alert",
            strong("Need ≥2 numeric phenotypes. "),
            sprintf(
              "%d plant row(s) after sidebar (%d before), but only %d eligible trait column(s).",
              n_after, nrow(ps_raw), length(eligible)
            )
          )
        } else {
          tags$div(
            class = "alert alert-info",
            role = "alert",
            sprintf(
              "%d plant rows after sidebar filters (%d before), %d eligible traits for PCA.",
              n_after, nrow(ps_raw), length(eligible)
            )
          )
        }
      }
    })

    observe({
      ps <- plants_sidebar()
      if (is.null(ps) || nrow(ps) < 1L) {
        ps <- .sf_dr_read_plants(plant_summaries)
      }
      choices <- .sf_dr_eligible_traits(ps)
      if (length(choices) < 1L) {
        updateSelectInput(session, "dr_traits", choices = character(0))
        return()
      }
      updateSelectInput(session, "dr_traits",
                        choices = choices,
                        selected = choices)
    })

    # Last Calculate PCA result (success or soft-fail); cleared when Growth Summaries changes
    dr_cache <- reactiveVal(NULL)

    observeEvent(
      plant_summaries(),
      {
        dr_cache(NULL)
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    .sf_dr_idle_msg <- function(what = "PCA results") {
      DT::datatable(
        data.frame(Message = paste0("Click Calculate PCA to show ", what, ".")),
        rownames = FALSE,
        options = list(dom = "t")
      )
    }

    .sf_dr_soft_fail <- function(reason) {
      DT::datatable(
        data.frame(Message = reason %||% "PCA failed"),
        rownames = FALSE,
        options = list(dom = "t")
      )
    }

    # Heavy PCA only on button; plant rows filtered by left sidebar
    dr_result <- eventReactive(input$calculate_dr, {
      ps <- plants_sidebar()
      # shiny:: — jsonlite::validate is attached via tidyverse and masks shiny::validate
      shiny::validate(
        shiny::need(
          !is.null(ps) && nrow(ps) > 0L,
          "No plants after left-sidebar filters. Widen cultivar/treatment or re-run Growth Summaries."
        ),
        shiny::need(
          length(input$dr_traits) >= 2L,
          "Select at least 2 phenotype columns."
        )
      )

      n_comp <- input$n_comp
      if (is.null(n_comp) || length(n_comp) < 1L || is.na(n_comp)) {
        n_comp <- NULL
      } else {
        n_comp <- as.integer(n_comp)
      }

      sidebar_snap <- .sf_sidebar_plant_snapshot(combined_inputs)

      withProgress(message = "Running plant × trait PCA…", value = 0.2, {
        incProgress(0.3, detail = paste(length(input$dr_traits), "traits | n=", nrow(ps)))
        out <- reduceTraitsPCA(
          plant_df = ps,
          trait_cols = input$dr_traits,
          n_comp = n_comp,
          scale = isTRUE(input$scale_traits),
          # GS multi-trait panel must be trait__auc / __peak / … (not raw scans)
          require_named_temporal_metrics = TRUE
        )
        if (is.list(out) && is.list(out$meta)) {
          out$meta$sidebar <- sidebar_snap
        }
        # Keep failure results so panels show soft-fail reason (not idle copy)
        dr_cache(out)
        incProgress(0.5, detail = "Done")
        out
      })
    })

    # Export PCs only on successful, still-aligned cache (failures stay in UI cache only)
    plant_with_pcs <- reactive({
      res <- dr_cache()
      if (is.null(res) || isFALSE(res$meta$success)) {
        return(NULL)
      }
      ps <- .sf_dr_read_plants(plant_summaries)
      pcs <- res$plant_df_with_pcs
      if (!.sf_dr_pcs_aligned(ps, pcs)) {
        return(NULL)
      }
      pcs
    })

    .sf_dr_fail_reason <- function(res) {
      if (!is.null(res$meta$reason) && nzchar(as.character(res$meta$reason))) {
        as.character(res$meta$reason)
      } else {
        "PCA failed"
      }
    }

    output$scree_plot <- renderPlot({
      res <- dr_cache()
      if (is.null(res)) {
        return(statfarmerEmptyPlot("Click Calculate PCA to draw the scree plot."))
      }
      if (isFALSE(res$meta$success)) {
        return(statfarmerEmptyPlot(.sf_dr_fail_reason(res)))
      }
      createDrScreePlot(res$variance)
    }, bg = "white")

    output$scores_plot <- renderPlot({
      res <- dr_cache()
      if (is.null(res)) {
        return(statfarmerEmptyPlot("Click Calculate PCA to draw plant PC scores."))
      }
      if (isFALSE(res$meta$success)) {
        return(statfarmerEmptyPlot(.sf_dr_fail_reason(res)))
      }
      color_by <- input$dr_color_by
      if (is.null(color_by) || identical(color_by, "none")) color_by <- NULL
      createDrScoresPlot(
        res$plant_df_with_pcs,
        color_by = color_by,
        title = "Plant PC scores (PC1 vs PC2)"
      )
    }, bg = "white")

    output$loadings_plot <- renderPlot({
      res <- dr_cache()
      if (is.null(res)) {
        return(statfarmerEmptyPlot("Click Calculate PCA to draw trait loadings."))
      }
      if (isFALSE(res$meta$success)) {
        return(statfarmerEmptyPlot(.sf_dr_fail_reason(res)))
      }
      createDrLoadingsPlot(res$loadings)
    }, bg = "white")

    .sf_dr_save_svg <- function(file, plot_obj) {
      ggplot2::ggsave(
        file, plot = plot_obj,
        width = 180, height = 112, units = "mm", bg = "white",
        device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
      )
    }

    output$scree_svg <- downloadHandler(
      filename = function() paste0("trait_dr_scree_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        .sf_dr_save_svg(file, createDrScreePlot(res$variance))
      }
    )
    output$scores_svg <- downloadHandler(
      filename = function() paste0("trait_dr_scores_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        color_by <- input$dr_color_by
        if (is.null(color_by) || identical(color_by, "none")) color_by <- NULL
        .sf_dr_save_svg(
          file,
          createDrScoresPlot(
            res$plant_df_with_pcs,
            color_by = color_by,
            title = "Plant PC scores (PC1 vs PC2)"
          )
        )
      }
    )
    output$loadings_svg <- downloadHandler(
      filename = function() paste0("trait_dr_loadings_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        .sf_dr_save_svg(file, createDrLoadingsPlot(res$loadings))
      }
    )

    output$variance_table <- DT::renderDataTable({
      res <- dr_cache()
      if (is.null(res)) {
        return(.sf_dr_idle_msg("variance table"))
      }
      if (isFALSE(res$meta$success)) {
        return(.sf_dr_soft_fail(res$meta$reason))
      }
      DT::datatable(
        res$variance %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$loadings_table <- DT::renderDataTable({
      res <- dr_cache()
      if (is.null(res)) {
        return(.sf_dr_idle_msg("loadings"))
      }
      if (isFALSE(res$meta$success)) {
        return(.sf_dr_soft_fail(res$meta$reason))
      }
      DT::datatable(
        res$loadings %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$scores_table <- DT::renderDataTable({
      res <- dr_cache()
      if (is.null(res)) {
        return(.sf_dr_idle_msg("plant scores"))
      }
      if (isFALSE(res$meta$success)) {
        return(.sf_dr_soft_fail(res$meta$reason))
      }
      df <- res$plant_df_with_pcs
      id_cols <- intersect(c("unit", "V.T.R", "cultivar", "treatment"), names(df))
      pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
      show <- df[, c(id_cols, pc_cols), drop = FALSE]
      DT::datatable(
        show %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$cluster_table <- DT::renderDataTable({
      res <- dr_cache()
      if (is.null(res)) {
        return(.sf_dr_idle_msg("trait clusters"))
      }
      if (isFALSE(res$meta$success)) {
        return(.sf_dr_soft_fail(res$meta$reason))
      }
      DT::datatable(
        res$cluster,
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$loadings_download <- downloadHandler(
      filename = function() paste0("trait_dr_loadings_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        write.csv(res$loadings, file, row.names = FALSE)
      }
    )

    output$scores_download <- downloadHandler(
      filename = function() paste0("trait_dr_plant_scores_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        write.csv(res$plant_df_with_pcs, file, row.names = FALSE)
      }
    )

    output$cluster_download <- downloadHandler(
      filename = function() paste0("trait_dr_clusters_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- dr_cache()
        req(res, isTRUE(res$meta$success))
        write.csv(res$cluster, file, row.names = FALSE)
      }
    )

    list(plant_with_pcs = plant_with_pcs, dr_result = dr_result)
  })
}
