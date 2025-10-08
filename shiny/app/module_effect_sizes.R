# Effect Sizes Module
# Calculates Cohen's d, partial eta-squared, and other effect size measures

#' Effect Sizes UI
#' @param id Module ID
effectSizesUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Effect Sizes",
    fluidRow(
      column(12,
        h3("Effect Size Estimates"),
        p("Quantify the magnitude of differences between groups, independent of sample size."),
        hr()
      )
    ),
    
    fluidRow(
      column(4,
        wellPanel(
          h4("Settings"),
          selectInput(ns("es_factor"), "Factor for effect size:", 
                     choices = c(), multiple = FALSE),
          checkboxInput(ns("es_cohens_d"), "Cohen's d (pairwise)", value = TRUE),
          checkboxInput(ns("es_partial_eta"), "Partial eta-squared (ANOVA)", value = TRUE),
          checkboxInput(ns("es_omega"), "Omega-squared", value = FALSE),
          hr(),
          actionButton(ns("calculate_es"), "Calculate Effect Sizes", 
                      class = "btn-primary", icon = icon("calculator"))
        )
      ),
      
      column(8,
        tabsetPanel(
          tabPanel("Cohen's d",
            br(),
            p("Cohen's d measures standardized difference between two means."),
            p("Rules of thumb: |d| = 0.2 (small), 0.5 (medium), 0.8 (large)"),
            DT::dataTableOutput(ns("cohens_d_table")),
            br(),
            downloadButton(ns("cohens_d_download"), "Download Cohen's d (CSV)")
          ),
          tabPanel("Eta-squared",
            br(),
            p("Partial eta-squared (η²) measures proportion of variance explained by a factor."),
            p("Rules of thumb: η² = 0.01 (small), 0.06 (medium), 0.14 (large)"),
            DT::dataTableOutput(ns("eta_squared_table")),
            br(),
            downloadButton(ns("eta_squared_download"), "Download Eta-squared (CSV)")
          ),
          tabPanel("Interpretation",
            br(),
            h4("Effect Size Interpretation Guide"),
            tags$ul(
              tags$li(strong("Cohen's d:"), " Standardized mean difference. Independent of sample size."),
              tags$li(strong("Partial η²:"), " Variance explained by factor after controlling for other factors."),
              tags$li(strong("Omega-squared (ω²):"), " Less biased estimate than eta-squared."),
              tags$li(strong("Confidence intervals:"), " Show uncertainty in effect size estimates.")
            ),
            br(),
            h5("When to use:"),
            tags$ul(
              tags$li("Small p-values with small effects → statistically but not practically significant"),
              tags$li("Large effects with large p-values → may need more data"),
              tags$li("Comparing across studies → use standardized effect sizes")
            )
          )
        )
      )
    )
  )
}

#' Effect Sizes Server
#' @param id Module ID
#' @param filteredData Reactive filtered data
#' @param analysisResults Reactive analysis results
#' @param combined_inputs Reactive combined inputs
effectSizesServer <- function(id, filteredData, analysisResults, combined_inputs) {
  moduleServer(id, function(input, output, session) {
    
    # Update factor choices
    observe({
      req(combined_inputs$anova_factors)
      updateSelectInput(session, "es_factor", 
                       choices = combined_inputs$anova_factors,
                       selected = if ("treatment" %in% combined_inputs$anova_factors) "treatment" else combined_inputs$anova_factors[1])
    })
    
    # Calculate effect sizes
    effect_sizes <- eventReactive(input$calculate_es, {
      req(filteredData(), analysisResults(), input$es_factor)
      
      data <- filteredData()
      factor <- input$es_factor
      outcome <- combined_inputs$out_variables
      
      cat("DEBUG effectSizes: factor =", factor, ", outcome =", outcome, "\n")
      cat("DEBUG effectSizes: data dims =", nrow(data), "x", ncol(data), "\n")
      cat("DEBUG effectSizes: factor in data?", factor %in% names(data), "\n")
      cat("DEBUG effectSizes: outcome in data?", outcome %in% names(data), "\n")
      
      results <- list()
      
      # Cohen's d for all pairwise comparisons
      if (input$es_cohens_d) {
        tryCatch({
          levels <- unique(data[[factor]])
          cat("DEBUG effectSizes: factor levels =", paste(levels, collapse=", "), "\n")
          comparisons <- combn(levels, 2, simplify = FALSE)
          cat("DEBUG effectSizes: n comparisons =", length(comparisons), "\n")
          
          cohens_d_list <- lapply(comparisons, function(pair) {
            # Extract as vector, not data frame
            group1 <- data[data[[factor]] == pair[1], outcome, drop = TRUE]
            group2 <- data[data[[factor]] == pair[2], outcome, drop = TRUE]
            
            m1 <- mean(group1, na.rm = TRUE)
            m2 <- mean(group2, na.rm = TRUE)
            sd1 <- sd(group1, na.rm = TRUE)
            sd2 <- sd(group2, na.rm = TRUE)
            n1 <- sum(!is.na(group1))
            n2 <- sum(!is.na(group2))
            
            # Pooled SD
            pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
            d <- (m1 - m2) / pooled_sd
            
            # 95% CI for d (approximate)
            se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
            ci_lower <- d - 1.96 * se_d
            ci_upper <- d + 1.96 * se_d
            
            # Interpretation
            interp <- ifelse(abs(d) < 0.2, "negligible",
                           ifelse(abs(d) < 0.5, "small",
                                 ifelse(abs(d) < 0.8, "medium", "large")))
            
            data.frame(
              comparison = paste(pair[1], "vs", pair[2]),
              cohens_d = d,
              ci_lower = ci_lower,
              ci_upper = ci_upper,
              interpretation = interp,
              stringsAsFactors = FALSE
            )
          })
          
          results$cohens_d <- do.call(rbind, cohens_d_list)
          cat("DEBUG effectSizes: Cohen's d calculated, rows =", nrow(results$cohens_d), "\n")
        }, error = function(e) {
          cat("DEBUG effectSizes: Cohen's d ERROR:", e$message, "\n")
          results$cohens_d <- data.frame(error = paste("Error calculating Cohen's d:", e$message))
        })
      }
      
      # Partial eta-squared from ANOVA
      if (input$es_partial_eta && !is.null(analysisResults()$anova$success) && analysisResults()$anova$success) {
        tryCatch({
          # For lmer models, partial eta-squared is not standard; report R^2 (marginal/conditional)
          model_type <- analysisResults()$anova$model_type
          if (!is.null(model_type) && identical(model_type, 'lmer')) {
            if (!requireNamespace('MuMIn', quietly = TRUE)) stop('MuMIn package is required for R-squared (GLMM)')
            r2 <- MuMIn::r.squaredGLMM(analysisResults()$anova$model)
            results$eta_squared <- data.frame(
              metric = c('R2_marginal', 'R2_conditional'),
              value = c(unname(r2[1]), unname(r2[2])),
              note = c('Variance explained by fixed effects', 'Variance explained by fixed + random effects'),
              stringsAsFactors = FALSE
            )
          } else {
            an <- analysisResults()$anova$anova_table
            # Robustly detect sums of squares and term columns (support broom::tidy names too)
            ss_col <- c('Sum Sq', 'Sum.Sq', 'SumSq', 'SS', 'sumsq')
            ss_col <- ss_col[ss_col %in% names(an)]
            # If not found, fallback to base anova() on the model
            if (length(ss_col) == 0 && !is.null(analysisResults()$anova$model)) {
              base_an <- try(as.data.frame(anova(analysisResults()$anova$model)), silent = TRUE)
              if (!inherits(base_an, 'try-error')) {
                an <- base_an
                ss_col <- c('Sum Sq', 'Sum.Sq', 'SumSq', 'SS')
                ss_col <- ss_col[ss_col %in% names(an)]
              }
            }
            if (length(ss_col) == 0) stop('No Sum of Squares column found in ANOVA table')
            ss_col <- ss_col[1]

            # Determine term names and residual row
            if (!is.null(rownames(an)) && any(rownames(an) == 'Residuals')) {
              terms_vec <- rownames(an)
              residual_ss <- an['Residuals', ss_col][[1]]
              term_rows <- terms_vec != 'Residuals'
              term_names <- terms_vec[term_rows]
              term_ss <- an[term_rows, ss_col][[1]]
            } else if ('term' %in% names(an)) {
              term_names <- an$term[an$term != 'Residuals']
              term_ss <- an[[ss_col]][an$term != 'Residuals']
              residual_ss <- an[[ss_col]][an$term == 'Residuals'][1]
            } else if ('Term' %in% names(an)) {
              term_names <- an$Term[an$Term != 'Residuals']
              term_ss <- an[[ss_col]][an$Term != 'Residuals']
              residual_ss <- an[[ss_col]][an$Term == 'Residuals'][1]
            } else {
              stop('Cannot identify term names in ANOVA table')
            }

            eta_vals <- term_ss / (term_ss + residual_ss)
            interpretation <- ifelse(eta_vals < 0.01, 'negligible',
                                      ifelse(eta_vals < 0.06, 'small',
                                             ifelse(eta_vals < 0.14, 'medium', 'large')))
            results$eta_squared <- data.frame(
              term = term_names,
              partial_eta_sq = eta_vals,
              interpretation = interpretation,
              stringsAsFactors = FALSE
            )
          }
        }, error = function(e) {
          results$eta_squared <- data.frame(error = paste("Error calculating eta-squared:", e$message))
        })
      }
      
      cat("DEBUG effectSizes: Returning results. Cohen's d?", !is.null(results$cohens_d), 
          ", Eta?", !is.null(results$eta_squared), "\n")
      results
    })
    
    # Render Cohen's d table
    output$cohens_d_table <- DT::renderDataTable({
      req(effect_sizes())
      cat("DEBUG render cohens_d: effect_sizes() exists\n")
      if (!is.null(effect_sizes()$cohens_d)) {
        cat("DEBUG render cohens_d: has data, rows =", nrow(effect_sizes()$cohens_d), "\n")
        DT::datatable(
          effect_sizes()$cohens_d %>%
            mutate(across(tidyselect::where(is.numeric), ~signif(., 3))),
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        )
      } else {
        cat("DEBUG render cohens_d: NULL, showing message\n")
        DT::datatable(data.frame(Message = "Click 'Calculate Effect Sizes' to compute"))
      }
    })
    
    # Render eta-squared table
    output$eta_squared_table <- DT::renderDataTable({
      req(effect_sizes())
      if (!is.null(effect_sizes()$eta_squared)) {
        DT::datatable(
          effect_sizes()$eta_squared %>%
            mutate(across(tidyselect::where(is.numeric), ~signif(., 3))),
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        )
      } else {
        DT::datatable(data.frame(Message = "Click 'Calculate Effect Sizes' to compute"))
      }
    })

    # Downloads
    output$cohens_d_download <- downloadHandler(
      filename = function() paste0("cohens_d_", Sys.Date(), ".csv"),
      content = function(file) {
        req(effect_sizes(), effect_sizes()$cohens_d)
        write.csv(effect_sizes()$cohens_d, file, row.names = FALSE)
      }
    )

    output$eta_squared_download <- downloadHandler(
      filename = function() paste0("eta_squared_", Sys.Date(), ".csv"),
      content = function(file) {
        req(effect_sizes(), effect_sizes()$eta_squared)
        write.csv(effect_sizes()$eta_squared, file, row.names = FALSE)
      }
    )
    
  })
}

