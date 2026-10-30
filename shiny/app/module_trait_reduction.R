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
        p("PCA on plant × selected traits (never genotype rows or design factors). PC scores attach to the same plant rows for BLUP."),
        uiOutput(ns("dr_guidance")),
        hr()
      )
    ),

    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          selectInput(ns("dr_traits"), "Traits (plant phenotypes):",
                      choices = c(), multiple = TRUE),
          numericInput(ns("n_comp"), "Number of PCs (blank = all):",
                       value = NA, min = 1, max = 50, step = 1),
          checkboxInput(ns("scale_traits"), "Scale traits before PCA", value = TRUE),
          helpText(
            tags$small(
              tags$strong("Rows:"), " plants from Growth Summaries (one row per plant).", tags$br(),
              tags$strong("Columns:"), " phenotypes only — cultivar/treatment stay out of PCA.", tags$br(),
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
          tabPanel("Scree",
            br(),
            p("Proportion of variance per PC. Click Calculate PCA to fill this panel."),
            plotOutput(ns("scree_plot"), height = "360px"),
            br(),
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
              tags$li(strong("Rows:"), " plants only (Growth Summaries collapse)."),
              tags$li(strong("Columns:"), " user-selected numeric phenotypes — not cultivar/treatment."),
              tags$li(strong("Output:"), " PC1… on the same plant table for the BLUP trait chooser."),
              tags$li(strong("Refresh:"), " re-running Growth Summaries clears cached PCs — Calculate PCA again."),
              tags$li(strong("Anti-goal:"), " never mass-BLUP-then-PCA; never genotype×trait mean matrix."),
              tags$li(strong("Parked:"), " multi-PlantEye wide plant×traits frame (beyond GS curve metrics) is follow-up.")
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

    output$dr_guidance <- renderUI({
      ps <- .sf_dr_read_plants(plant_summaries)
      if (is.null(ps)) {
        tags$div(
          class = "alert alert-warning",
          role = "alert",
          strong("Growth Summaries required. "),
          "Open Growth Summaries and click Calculate Summaries so Trait DR has one row per plant."
        )
      } else {
        eligible <- .sf_dr_eligible_traits(ps)
        if (length(eligible) < 2L) {
          tags$div(
            class = "alert alert-warning",
            role = "alert",
            strong("Need ≥2 numeric phenotypes. "),
            sprintf(
              "%d plant row(s) ready, but only %d eligible trait column(s). Enable more Growth Summaries metrics (mean/AUC/peak/rate) or provide a wider plant-level table.",
              nrow(ps), length(eligible)
            )
          )
        } else {
          tags$div(
            class = "alert alert-info",
            role = "alert",
            sprintf(
              "%d plant rows, %d eligible traits ready for plant×trait PCA.",
              nrow(ps), length(eligible)
            )
          )
        }
      }
    })

    observe({
      ps <- .sf_dr_read_plants(plant_summaries)
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

    # Heavy PCA only on button (shiny-r-optimize)
    dr_result <- eventReactive(input$calculate_dr, {
      ps <- .sf_dr_read_plants(plant_summaries)
      validate(
        need(
          !is.null(ps),
          "Run Growth Summaries first (Calculate Summaries) before Trait DR."
        ),
        need(
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

      withProgress(message = "Running plant × trait PCA…", value = 0.2, {
        incProgress(0.3, detail = paste(length(input$dr_traits), "traits"))
        out <- reduceTraitsPCA(
          plant_df = ps,
          trait_cols = input$dr_traits,
          n_comp = n_comp,
          scale = isTRUE(input$scale_traits)
        )
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

    output$scree_plot <- renderPlot({
      res <- dr_cache()
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Click Calculate PCA to draw the scree plot.", cex = 1.1)
        return(invisible())
      }
      if (isFALSE(res$meta$success)) {
        plot.new()
        text(0.5, 0.5, res$meta$reason %||% "PCA failed", cex = 1.2)
        return(invisible())
      }
      v <- res$variance
      req(v, nrow(v) > 0L)
      op <- par(mar = c(5, 4, 2, 1))
      on.exit(par(op), add = TRUE)
      barplot(
        v$prop_var,
        names.arg = v$component,
        ylab = "Proportion of variance",
        xlab = "Component",
        main = "Scree (plant × traits PCA)",
        col = "#3498db",
        border = NA
      )
    })

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
