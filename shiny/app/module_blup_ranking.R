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
          tabPanel("Answer",
            br(),
            uiOutput(ns("blup_answer_panel")),
            hr(),
            p("Follow FarmaClip (📎 bottom-right): Growth Summaries → Trait DR → BLUP on PC1/PC2 with a check cultivar.")
          ),
          tabPanel("Curves & Tukey",
            br(),
            p("Curves use the full time series. Violins + Tukey use ONE time slice only — never a raw mean across dates."),
            selectInput(ns("followup_trait"), "Follow-up trait:", choices = c(), multiple = FALSE),
            selectInput(ns("followup_time"), "Violin/Tukey time slice:", choices = c(), multiple = FALSE),
            helpText(tags$small("Default = last cluster. Change slice to compare genotypes at that developmental state.")),
            actionButton(ns("run_followup"), "Plot curves + Tukey letters",
                         class = "btn-primary", icon = icon("chart-area")),
            hr(),
            uiOutput(ns("followup_slice_note")),
            plotOutput(ns("followup_curve"), height = "420px"),
            plotOutput(ns("followup_violin"), height = "420px"),
            br(),
            DT::dataTableOutput(ns("followup_tukey_letters"))
          ),
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
#' @param filteredData Optional reactive observation-level data (curves / Tukey follow-up)
#' @param dr_result Optional reactive Trait DR result (suggests follow-up traits from loadings)
blupRankingServer <- function(id, plant_summaries, combined_inputs,
                              filteredData = NULL, dr_result = NULL) {
  moduleServer(id, function(input, output, session) {

    # Plant table after left-sidebar cultivar/treatment filters (same as main plot)
    plants_sidebar <- reactive({
      ps <- .sf_blup_read_plants(plant_summaries)
      .sf_filter_plants_by_sidebar(ps, combined_inputs)
    })

    output$blup_guidance <- renderUI({
      ps_raw <- .sf_blup_read_plants(plant_summaries)
      if (is.null(ps_raw)) {
        tags$div(
          class = "alert alert-warning",
          role = "alert",
          strong("Growth Summaries required. "),
          "Open the Growth Summaries tab and click Calculate Summaries to collapse plants over the selected time clusters. BLUP will not run on multi-time raw rows."
        )
      } else {
        ps <- plants_sidebar()
        n_traits <- length(.sf_blup_trait_choices(ps_raw))
        n_after <- if (is.null(ps)) 0L else nrow(ps)
        tags$div(
          class = "alert alert-info",
          role = "alert",
          sprintf(
            "%d plant rows after sidebar filters (%d before; %d selectable traits). Cultivar/treatment from the left bar apply to fit + plots. Click Calculate BLUP.",
            n_after, nrow(ps_raw), n_traits
          )
        )
      }
    })

    observe({
      ps <- plants_sidebar()
      if (is.null(ps) || nrow(ps) < 1L) {
        ps <- .sf_blup_read_plants(plant_summaries)
      }
      choices <- .sf_blup_trait_choices(ps)
      if (length(choices) < 1L) {
        updateSelectInput(session, "blup_trait", choices = character(0))
        return()
      }
      # Prefer named GS metrics (multi-trait panel), then legacy single-trait names, then PCs
      preferred <- c(
        grep("__auc$", choices, value = TRUE),
        grep("__last$", choices, value = TRUE),
        grep("__peak$", choices, value = TRUE),
        intersect(c("auc", "last", "peak", "growth_rate", "mean"), choices),
        grep("^PC[0-9]+$", choices, value = TRUE)
      )
      preferred <- unique(preferred)
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

    # Heavy fit only on button (shiny-r-optimize); plants filtered by left sidebar
    blup_result <- eventReactive(input$calculate_blup, {
      ps <- plants_sidebar()
      # shiny:: — jsonlite::validate is attached via tidyverse and masks shiny::validate
      shiny::validate(
        shiny::need(
          !is.null(ps) && nrow(ps) > 0L,
          "No plants after left-sidebar filters. Widen cultivar/treatment selection or re-run Growth Summaries."
        ),
        shiny::need(
          shiny::isTruthy(input$blup_trait) && input$blup_trait %in% names(ps),
          "Select a plant-level trait column."
        ),
        shiny::need(
          "cultivar" %in% names(ps),
          "Plant summaries missing cultivar column."
        )
      )

      checks <- .sf_blup_parse_checks(input$checks)
      # Checks only among genotypes still present after sidebar filter
      if (!is.null(checks) && "cultivar" %in% names(ps)) {
        checks <- intersect(checks, unique(as.character(ps$cultivar)))
        if (length(checks) < 1L) checks <- NULL
      }
      top_k <- as.integer(input$top_k)
      if (is.na(top_k) || top_k < 1L) {
        top_k <- 20L
      }

      sidebar_snap <- .sf_sidebar_plant_snapshot(combined_inputs)

      withProgress(message = "Fitting genotype BLUP…", value = 0.2, {
        incProgress(0.3, detail = paste("Trait:", input$blup_trait, "| n=", nrow(ps)))
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
        if (is.list(out) && is.list(out$meta)) {
          out$meta$sidebar <- sidebar_snap
          out$meta$n_plants_sidebar <- nrow(ps)
        }
        incProgress(0.5, detail = "Done")
        out
      })
    })

    .sf_blup_sidebar_stale <- function(res) {
      if (is.null(res) || is.null(res$meta$sidebar)) return(FALSE)
      cur <- .sf_sidebar_plant_snapshot(combined_inputs)
      !identical(cur$cultivars, res$meta$sidebar$cultivars) ||
        !identical(cur$treatments, res$meta$sidebar$treatments)
    }

    .sf_blup_plot_geno <- function(res) {
      .sf_blup_filter_geno_table(res$genotype_table, combined_inputs)
    }
    .sf_blup_plot_trt <- function(res) {
      .sf_blup_filter_trt_table(res$treatment_table, combined_inputs)
    }
    .sf_blup_plot_dunnett <- function(res) {
      .sf_blup_filter_dunnett(res$meta$dunnett, combined_inputs)
    }

    .sf_blup_fail_msg <- function(res) {
      if (!is.null(res$reason) && nzchar(as.character(res$reason))) {
        as.character(res$reason)
      } else {
        "BLUP failed"
      }
    }

    output$blup_answer_panel <- renderUI({
      if (!isTruthy(input$calculate_blup) || input$calculate_blup < 1L) {
        return(tags$div(
          class = "alert alert-info",
          "Calculate BLUP (preferably on PC1 after Trait DR) to fill this answer panel."
        ))
      }
      res <- blup_result()
      req(res)
      if (isFALSE(res$success)) {
        return(tags$div(class = "alert alert-warning", .sf_blup_fail_msg(res)))
      }
      sig <- .sf_dunnett_significant(.sf_blup_plot_dunnett(res), alpha = 0.05)
      geno <- .sf_blup_plot_geno(res)
      trait <- res$meta$trait %||% "trait"
      checks <- res$meta$checks
      tagList(
        tags$h4(sprintf("Results for %s", trait)),
        if (!is.null(checks) && length(checks)) {
          tags$p(em(sprintf("Dunnett check(s): %s", paste(checks, collapse = ", "))))
        },
        if (is.data.frame(sig) && nrow(sig) > 0L) {
          tags$div(
            class = "alert alert-success",
            tags$strong(sprintf("%d significant contrast(s) (FDR p < 0.05)", nrow(sig))),
            DT::datatable(
              sig %>% dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
              options = list(dom = "t", pageLength = 12),
              rownames = FALSE
            )
          )
        } else {
          tags$div(
            class = "alert alert-secondary",
            "No Dunnett FDR hits at α=0.05. Top BLUP ranks (below) are still the ranking answer."
          )
        },
        tags$h5("Top genotype BLUPs"),
        DT::datatable(
          utils::head(geno, 12L) %>%
            dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 4))),
          options = list(dom = "t", pageLength = 12),
          rownames = FALSE
        )
      )
    })

    .sf_blup_rank_title <- function(res) {
      trait <- res$meta$trait
      base <- paste0(
        "Genotype BLUP ranks",
        if (!is.null(trait) && nzchar(as.character(trait))) paste0(" — ", trait) else ""
      )
      if (.sf_blup_sidebar_stale(res)) {
        paste0(base, " (sidebar changed — re-run Calculate for refit)")
      } else {
        paste0(base, " (sidebar filters)")
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
      createBlupRankPlot(
        .sf_blup_plot_geno(res),
        checks = res$meta$checks,
        title = .sf_blup_rank_title(res)
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
      createBlupDunnettPlot(
        .sf_blup_plot_dunnett(res),
        title = if (.sf_blup_sidebar_stale(res)) {
          "Dunnett vs check (FDR) — sidebar changed; re-run Calculate"
        } else {
          "Dunnett vs check (FDR) — sidebar cultivars"
        }
      )
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
      createBlupTreatmentPlot(
        .sf_blup_plot_trt(res),
        title = if (.sf_blup_sidebar_stale(res)) {
          "Treatment effects — sidebar changed; re-run Calculate"
        } else {
          "Treatment effects — sidebar treatments"
        }
      )
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
      createBlupRankPlot(
        .sf_blup_plot_geno(res),
        checks = res$meta$checks,
        title = .sf_blup_rank_title(res)
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
        .sf_blup_save_svg(file, createBlupDunnettPlot(.sf_blup_plot_dunnett(res)))
      }
    )
    output$blup_treatment_svg <- downloadHandler(
      filename = function() paste0("blup_treatment_", Sys.Date(), ".svg"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success))
        .sf_blup_save_svg(file, createBlupTreatmentPlot(.sf_blup_plot_trt(res)))
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
      d <- .sf_blup_plot_dunnett(res)
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
          data.frame(Message = "Dunnett table empty after sidebar cultivar filter."),
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
        d <- .sf_blup_plot_dunnett(res)
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
      tbl <- .sf_blup_plot_geno(res)
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
      tbl <- .sf_blup_plot_trt(res)
      if (is.null(tbl) || nrow(tbl) < 1L) {
        return(DT::datatable(
          data.frame(
            Message = "No treatment rows — see note above (1-level treatment or sidebar filter)."
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
        req(res, isTRUE(res$success))
        write.csv(.sf_blup_plot_geno(res), file, row.names = FALSE)
      }
    )

    output$treatment_download <- downloadHandler(
      filename = function() paste0("blup_treatment_", Sys.Date(), ".csv"),
      content = function(file) {
        res <- blup_result()
        req(res, isTRUE(res$success))
        write.csv(.sf_blup_plot_trt(res), file, row.names = FALSE)
      }
    )

    # ----- Curves & Tukey follow-up (obs-level; violin/Tukey = single time only) -----
    observe({
      fd <- tryCatch({
        if (is.null(filteredData)) NULL else filteredData()
      }, error = function(e) NULL)
      choices <- character(0)
      dr <- tryCatch({
        if (is.null(dr_result)) NULL else dr_result()
      }, error = function(e) NULL)
      if (!is.null(dr) && isTRUE(dr$meta$success)) {
        choices <- suggestFollowupTraits(dr$loadings, max_n = 8L)
      }
      if (!is.null(fd) && is.data.frame(fd)) {
        num <- names(fd)[vapply(fd, is.numeric, logical(1))]
        design <- c(
          "unit", "dbscan_cluster", "pos_x", "pos_y", "table_no", "repeat",
          "outlier", "time_numeric"
        )
        pool <- setdiff(num, design)
        if (length(choices)) {
          choices <- intersect(choices, pool)
          if (length(choices) < 1L) choices <- pool
        } else {
          choices <- pool
        }
        pref <- grep("biomass|ndvi|height|leaf_area|growth", choices, ignore.case = TRUE, value = TRUE)
        choices <- unique(c(pref, choices))
        tlev <- .sf_followup_time_levels(fd)
        if (length(tlev) > 0L) {
          updateSelectInput(
            session, "followup_time",
            choices = tlev,
            selected = tlev[[length(tlev)]]
          )
        } else {
          updateSelectInput(session, "followup_time", choices = character(0))
        }
      }
      if (length(choices) < 1L) {
        updateSelectInput(session, "followup_trait", choices = character(0))
      } else {
        updateSelectInput(session, "followup_trait", choices = choices, selected = choices[[1L]])
      }
    })

    followup_bundle <- eventReactive(input$run_followup, {
      fd <- tryCatch({
        if (is.null(filteredData)) NULL else filteredData()
      }, error = function(e) NULL)
      shiny::validate(
        shiny::need(!is.null(fd) && is.data.frame(fd) && nrow(fd) > 0L,
             "No filtered observation data — pick a processed project and filters in the left bar."),
        shiny::need(shiny::isTruthy(input$followup_trait) && input$followup_trait %in% names(fd),
             "Select a follow-up trait present in the filtered data.")
      )
      fd2 <- .sf_filter_plants_by_sidebar(fd, combined_inputs)
      if (is.null(fd2) || nrow(fd2) < 1L) fd2 <- fd
      tv <- input$followup_time
      if (!isTruthy(tv)) tv <- NULL
      buildTraitFollowupPlots(
        fd2, input$followup_trait,
        group = "cultivar",
        time_value = tv,
        run_tukey = TRUE
      )
    })

    output$followup_slice_note <- renderUI({
      if (!isTruthy(input$run_followup) || input$run_followup < 1L) return(NULL)
      b <- followup_bundle()
      req(b)
      tags$div(
        class = "alert alert-info",
        role = "alert",
        tags$strong("Time discipline: "),
        b$slice_note %||% "Single-time violin/Tukey.",
        " Growth curve above still shows the full temporal trajectory."
      )
    })

    output$followup_curve <- renderPlot({
      if (!isTruthy(input$run_followup) || input$run_followup < 1L) {
        return(statfarmerEmptyPlot("Pick a trait + time slice, then Plot curves + Tukey letters."))
      }
      b <- followup_bundle()
      req(b)
      b$curve
    }, bg = "white")

    output$followup_violin <- renderPlot({
      if (!isTruthy(input$run_followup) || input$run_followup < 1L) {
        return(statfarmerEmptyPlot("Single-time Tukey violins appear after Plot curves + Tukey letters."))
      }
      b <- followup_bundle()
      req(b)
      b$violin
    }, bg = "white")

    output$followup_tukey_letters <- DT::renderDataTable({
      if (!isTruthy(input$run_followup) || input$run_followup < 1L) {
        return(DT::datatable(
          data.frame(Message = "Tukey CLD letters (single time slice) fill after follow-up plot."),
          rownames = FALSE, options = list(dom = "t")
        ))
      }
      b <- followup_bundle()
      req(b)
      if (is.null(b$tukey) || !isTRUE(b$tukey$success) || is.null(b$tukey$letters)) {
        msg <- if (!is.null(b$tukey$reason)) b$tukey$reason else "Tukey letters unavailable"
        return(DT::datatable(data.frame(Message = msg), rownames = FALSE, options = list(dom = "t")))
      }
      DT::datatable(
        b$tukey$letters,
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    invisible(list(blup_result = blup_result, followup_bundle = followup_bundle))
  })
}
