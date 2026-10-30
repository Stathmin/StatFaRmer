# BLUP Ranking Module
# Genotype BLUP on plant-collapsed traits (Growth Summaries → fitBlupModels)

.sf_blup_design_cols <- c(
  "unit", "V.T.R", "cultivar", "treatment", "repeat",
  "T:X:Y", "table_no", "pos_x", "pos_y",
  "dbscan_cluster", "timestamp_group", "n_obs", "n_times"
)

#' Safe read of plant_summaries reactive (NULL if not calculated yet).
.sf_blup_read_plants <- function(plant_summaries) {
  tryCatch({
    ps <- plant_summaries()
    if (is.null(ps) || !is.data.frame(ps) || nrow(ps) < 1L) {
      return(NULL)
    }
    ps
  }, error = function(e) NULL)
}

#' Numeric trait candidates from plant-level summaries (excludes design cols).
.sf_blup_trait_choices <- function(ps) {
  if (is.null(ps) || !is.data.frame(ps)) {
    return(character(0))
  }
  num <- names(ps)[vapply(ps, is.numeric, logical(1))]
  setdiff(num, .sf_blup_design_cols)
}

#' Parse optional comma-separated check genotype levels.
.sf_blup_parse_checks <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(as.character(txt)))) {
    return(NULL)
  }
  parts <- trimws(unlist(strsplit(as.character(txt), ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) < 1L) {
    return(NULL)
  }
  parts
}

#' BLUP UI
#' @param id Module ID
blupRankingUI <- function(id) {
  ns <- NS(id)

  tabPanel(
    "BLUP",
    fluidRow(
      column(12,
        h3("Genotype BLUP Ranking"),
        p("Rank genotypes with plant-level BLUPs after collapsing over selected times in Growth Summaries."),
        uiOutput(ns("blup_guidance")),
        hr()
      )
    ),

    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          selectInput(ns("blup_trait"), "Plant-level trait:",
                      choices = c(), multiple = FALSE),
          checkboxInput(ns("apply_spatial"), "Spatial residualize (per-table X/Y)", value = TRUE),
          numericInput(ns("top_k"), "Top-k genotypes:", value = 20, min = 1, max = 500, step = 1),
          textInput(ns("checks"), "Check genotypes (comma-separated, optional):",
                    value = "", placeholder = "e.g. CheckA, CheckB"),
          helpText(
            tags$small(
              tags$strong("Trait:"), " mean / AUC / peak / rate from Growth Summaries, or PC1… after Trait DR.", tags$br(),
              tags$strong("Checks:"), " leave empty for top-k ranks only; set names to run Dunnett/FDR vs checks.", tags$br(),
              tags$strong("Calculate:"), " fits only when you click — nothing runs on idle tab open."
            )
          ),
          hr(),
          actionButton(ns("calculate_blup"), "Calculate BLUP",
                       class = "btn-primary", icon = icon("calculator"))
        )
      ),

      column(8,
        tabsetPanel(
          tabPanel("Plots",
            br(),
            p("Shared StatFaRmer theme (same as main distPlot). Ranks, Dunnett posthoc, and treatment EMMeans."),
            plotOutput(ns("blup_rank_plot"), height = "420px"),
            downloadButton(ns("blup_rank_svg"), "Save ranks SVG"),
            hr(),
            plotOutput(ns("blup_dunnett_plot"), height = "360px"),
            downloadButton(ns("blup_dunnett_svg"), "Save Dunnett SVG"),
            hr(),
            plotOutput(ns("blup_treatment_plot"), height = "320px"),
            downloadButton(ns("blup_treatment_svg"), "Save treatment SVG")
          ),
          tabPanel("Genotype ranks",
            br(),
            p("Rank 1 = highest genotype BLUP. SE from conditional variance when available."),
            DT::dataTableOutput(ns("genotype_table")),
            br(),
            downloadButton(ns("genotype_download"), "Download genotype ranks (CSV)")
          ),
          tabPanel("Dunnett",
            br(),
            p("Fixed-genotype Dunnett vs checks (FDR). Empty until checks are set and Calculate BLUP succeeds."),
            DT::dataTableOutput(ns("dunnett_table")),
            br(),
            downloadButton(ns("dunnett_download"), "Download Dunnett (CSV)")
          ),
          tabPanel("Treatment",
            br(),
            uiOutput(ns("treatment_message")),
            DT::dataTableOutput(ns("treatment_table")),
            br(),
            downloadButton(ns("treatment_download"), "Download treatment table (CSV)")
          ),
          tabPanel("Model",
            br(),
            h4("Formula"),
            verbatimTextOutput(ns("blup_formula")),
            h4("Variance components"),
            verbatimTextOutput(ns("blup_varcomp")),
            h4("Heritability (narrow-sense proxy)"),
            verbatimTextOutput(ns("blup_h2")),
            h4("Notes"),
            verbatimTextOutput(ns("blup_meta_notes"))
          ),
          tabPanel("About",
            br(),
            h4("BLUP sidecar"),
            tags$ul(
              tags$li(strong("Input:"), " one row per plant from Growth Summaries (selected times only)."),
              tags$li(strong("Spatial:"), " optional per-table X/Y residualize before lmer."),
              tags$li(strong("Model:"), " (1|genotype) + treatment FE (≥2 levels) + (1|repeat) + (1|table_no)."),
              tags$li(strong("Plots:"), " createBlupRankPlot / Dunnett / Treatment via shiny/common/plotting.R."),
              tags$li(strong("Single treatment:"), " genotype ranks still fill; Treatment panel stays empty until ≥2 levels."),
              tags$li(strong("Multi-time raw:"), " blocked — collapse in Growth Summaries first.")
            )
          )
        )
      )
    )
  )
}

#' BLUP Server
#' @param id Module ID
#' @param plant_summaries Reactive plant-level summaries from growthSummariesServer
#' @param combined_inputs Reactive combined inputs
blupRankingServer <- function(id, plant_summaries, combined_inputs) {
  moduleServer(id, function(input, output, session) {

    output$blup_guidance <- renderUI({
      ps <- .sf_blup_read_plants(plant_summaries)
      if (is.null(ps)) {
        tags$div(
          class = "alert alert-warning",
          role = "alert",
          strong("Growth Summaries required. "),
          "Open the Growth Summaries tab and click Calculate Summaries to collapse plants over the selected time clusters. BLUP will not run on multi-time raw rows."
        )
      } else {
        n_traits <- length(.sf_blup_trait_choices(ps))
        tags$div(
          class = "alert alert-info",
          role = "alert",
          sprintf(
            "%d plant rows ready from Growth Summaries (%d selectable trait columns). Click Calculate BLUP to fit.",
            nrow(ps), n_traits
          )
        )
      }
    })

    observe({
      ps <- .sf_blup_read_plants(plant_summaries)
      choices <- .sf_blup_trait_choices(ps)
      if (length(choices) < 1L) {
        updateSelectInput(session, "blup_trait", choices = character(0))
        return()
      }
      preferred <- intersect(c("mean", "auc", "peak", "growth_rate"), choices)
      selected <- if (length(preferred) > 0L) preferred[[1L]] else choices[[1L]]
      updateSelectInput(session, "blup_trait",
                        choices = choices,
                        selected = selected)
    })

    .sf_blup_idle_msg <- function() {
      DT::datatable(
        data.frame(Message = "Click Calculate BLUP to rank genotypes."),
        rownames = FALSE,
        options = list(dom = "t")
      )
    }

    # Heavy fit only on button (shiny-r-optimize)
    blup_result <- eventReactive(input$calculate_blup, {
      ps <- .sf_blup_read_plants(plant_summaries)
      validate(
        need(
          !is.null(ps),
          "Run Growth Summaries first (Calculate Summaries) before BLUP."
        ),
        need(
          isTruthy(input$blup_trait) && input$blup_trait %in% names(ps),
          "Select a plant-level trait column."
        ),
        need(
          "cultivar" %in% names(ps),
          "Plant summaries missing cultivar column."
        )
      )

      checks <- .sf_blup_parse_checks(input$checks)
      top_k <- as.integer(input$top_k)
      if (is.na(top_k) || top_k < 1L) {
        top_k <- 20L
      }

      withProgress(message = "Fitting genotype BLUP…", value = 0.2, {
        incProgress(0.3, detail = paste("Trait:", input$blup_trait))
        out <- fitBlupModels(
          data = ps,
          trait = input$blup_trait,
          genotype_var = "cultivar",
          treatment_var = "treatment",
          repeat_var = "repeat",
          table_var = "table_no",
          apply_spatial = isTRUE(input$apply_spatial),
          checks = checks,
          top_k = top_k
        )
        incProgress(0.5, detail = "Done")
        out
      })
    })

    .sf_blup_fail_msg <- function(res) {
      if (!is.null(res$reason) && nzchar(as.character(res$reason))) {
        as.character(res$reason)
      } else {
        "BLUP failed"
      }
    }

    output$blup_rank_plot <- renderPlot({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(statfarmerEmptyPlot("Click Calculate BLUP to draw genotype ranks."))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(statfarmerEmptyPlot(.sf_blup_fail_msg(res)))
      }
      trait <- res$meta$trait
      createBlupRankPlot(
        res$genotype_table,
        checks = res$meta$checks,
        title = paste0(
          "Genotype BLUP ranks",
          if (!is.null(trait) && nzchar(as.character(trait))) paste0(" — ", trait) else ""
        )
      )
    }, bg = "white")

    output$blup_dunnett_plot <- renderPlot({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(statfarmerEmptyPlot("Dunnett plot appears after Calculate BLUP (with checks)."))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(statfarmerEmptyPlot(.sf_blup_fail_msg(res)))
      }
      createBlupDunnettPlot(res$meta$dunnett)
    }, bg = "white")

    output$blup_treatment_plot <- renderPlot({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(statfarmerEmptyPlot("Treatment plot appears after Calculate BLUP."))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(statfarmerEmptyPlot(.sf_blup_fail_msg(res)))
      }
      createBlupTreatmentPlot(res$treatment_table)
    }, bg = "white")

    .sf_blup_save_svg <- function(file, plot_obj) {
      ggplot2::ggsave(
        file, plot = plot_obj,
        width = 180, height = 112, units = "mm", bg = "white",
        device = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
      )
    }

    .sf_blup_current_rank_plot <- function() {
      res <- blup_result()
      req(res, isTRUE(res$success))
      trait <- res$meta$trait
      createBlupRankPlot(
        res$genotype_table,
        checks = res$meta$checks,
        title = paste0(
          "Genotype BLUP ranks",
          if (!is.null(trait) && nzchar(as.character(trait))) paste0(" — ", trait) else ""
        )
      )
    }

    output$blup_rank_svg <- downloadHandler(
      filename = function() paste0("blup_ranks_", Sys.Date(), ".svg"),
      content = function(file) .sf_blup_save_svg(file, .sf_blup_current_rank_plot())
    )
    output$blup_dunnett_svg <- downloadHandler(
      filename = function() paste0("blup_dunnett_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success))
        .sf_blup_save_svg(file, createBlupDunnettPlot(res$meta$dunnett))
      }
    )
    output$blup_treatment_svg <- downloadHandler(
      filename = function() paste0("blup_treatment_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success))
        .sf_blup_save_svg(file, createBlupTreatmentPlot(res$treatment_table))
      }
    )

    output$dunnett_table <- DT::renderDataTable({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(DT::datatable(
          data.frame(Message = "Set check genotypes and click Calculate BLUP for Dunnett."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(DT::datatable(
          data.frame(Message = .sf_blup_fail_msg(res)),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      d <- res$meta$dunnett
      if (is.null(d)) {
        return(DT::datatable(
          data.frame(Message = "No Dunnett — leave checks empty means ranks only."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      if (is.list(d) && !is.data.frame(d)) {
        return(DT::datatable(
          data.frame(Message = d$message %||% "Dunnett unavailable"),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      if (!is.data.frame(d) || nrow(d) < 1L) {
        return(DT::datatable(
          data.frame(Message = "Dunnett table empty."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      DT::datatable(
        d %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$dunnett_download <- downloadHandler(
      filename = function() paste0("blup_dunnett_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success))
        d <- res$meta$dunnett
        req(is.data.frame(d), nrow(d) > 0L)
        write.csv(d, file, row.names = FALSE)
      }
    )

    output$genotype_table <- DT::renderDataTable({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(.sf_blup_idle_msg())
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(DT::datatable(
          data.frame(Message = .sf_blup_fail_msg(res)),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      tbl <- res$genotype_table
      req(tbl)
      DT::datatable(
        tbl %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$treatment_message <- renderUI({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(tags$p(em("Treatment estimates appear after Calculate BLUP (needs ≥2 treatment levels).")))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(tags$div(
          class = "alert alert-warning",
          role = "alert",
          res$reason %||% "BLUP failed"
        ))
      }
      msg <- res$meta$treatment_message
      tbl <- res$treatment_table
      empty_trt <- is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) < 1L
      if (empty_trt) {
        detail <- if (!is.null(msg) && nzchar(as.character(msg))) {
          as.character(msg)
        } else {
          "No treatment fixed-effect estimates for this fit."
        }
        return(tags$div(
          class = "alert alert-info",
          role = "alert",
          strong("Treatment panel empty. "),
          "Genotype ranks still appear under Genotype ranks when only one treatment level is present (e.g. soy). ",
          tags$span(detail)
        ))
      }
      if (!is.null(msg) && nzchar(as.character(msg))) {
        tags$p(em(as.character(msg)))
      } else {
        NULL
      }
    })

    output$treatment_table <- DT::renderDataTable({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(DT::datatable(
          data.frame(Message = "Click Calculate BLUP to estimate treatment effects."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(DT::datatable(
          data.frame(Message = res$reason %||% "BLUP failed"),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      tbl <- res$treatment_table
      if (is.null(tbl) || nrow(tbl) < 1L) {
        return(DT::datatable(
          data.frame(
            Message = "No treatment rows — see note above (typical when treatment has 1 level)."
          ),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      DT::datatable(
        tbl %>%
          dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$blup_formula <- renderText({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return("Click Calculate BLUP to show the fitted formula.")
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(res$reason %||% "BLUP failed")
      }
      as.character(res$formula %||% NA_character_)
    })

    output$blup_varcomp <- renderPrint({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        cat("Click Calculate BLUP to show variance components.")
        return(invisible())
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        cat(res$reason %||% "BLUP failed")
        return(invisible())
      }
      print(res$varcomp)
    })

    output$blup_h2 <- renderText({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return("Click Calculate BLUP to estimate h2.")
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(res$reason %||% "BLUP failed")
      }
      h2 <- res$h2
      if (is.null(h2) || length(h2) < 1L || is.na(h2)) {
        return("h2 not estimable")
      }
      sprintf("h2 = %.4f  (σ²_g / (σ²_g + σ²_e))", as.numeric(h2))
    })

    output$blup_meta_notes <- renderText({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return("Click Calculate BLUP to show fit notes.")
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(res$reason %||% "BLUP failed")
      }
      meta <- res$meta
      lines <- c(
        paste0("spatial_applied: ", isTRUE(meta$spatial_applied)),
        paste0("include_table_re: ", isTRUE(meta$include_table_re)),
        paste0("table_re_note: ", meta$table_re_note %||% NA_character_),
        paste0("singular: ", isTRUE(meta$singular)),
        paste0("n_genotypes: ", meta$n_genotypes %||% NA),
        paste0("n_obs: ", meta$n_obs %||% NA),
        paste0("y_col: ", meta$y_col %||% NA_character_)
      )
      if (!is.null(meta$dunnett_message) && nzchar(as.character(meta$dunnett_message))) {
        lines <- c(lines, paste0("dunnett: ", meta$dunnett_message))
      } else if (is.data.frame(meta$dunnett)) {
        lines <- c(lines, paste0("dunnett contrasts: ", nrow(meta$dunnett), " rows"))
      }
      paste(lines, collapse = "\n")
    })

    output$genotype_download <- downloadHandler(
      filename = function() paste0("blup_genotype_ranks_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success), res$genotype_table)
        write.csv(res$genotype_table, file, row.names = FALSE)
      }
    )

    output$treatment_download <- downloadHandler(
      filename = function() paste0("blup_treatment_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success), res$treatment_table)
        write.csv(res$treatment_table, file, row.names = FALSE)
      }
    )

    invisible(list(blup_result = blup_result))
  })
}
