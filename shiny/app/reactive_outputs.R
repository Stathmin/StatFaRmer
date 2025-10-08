# StatFaRmer Reactive Output Management
# All output rendering functions

#' Create reactive output functions
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param combined_inputs Reactive values for combined inputs
#' @param filteredData Reactive function for filtered data
#' @param projectData Reactive function for project data
#' @return List of reactive functions
createReactiveOutputs <- function(input, output, session, combined_inputs, filteredData, projectData) {
  
  # Display formula and model info
  output$formula <- renderUI({
    req(combined_inputs$out_variables)
    req(combined_inputs$anova_factors)
    req(filteredData())
    
    # Debug logging
    cat("DEBUG: out_variables =", combined_inputs$out_variables, "\n")
    cat("DEBUG: anova_factors =", paste(combined_inputs$anova_factors, collapse = ", "), "\n")
    
    # Analyze cardinality to preview model choice
    user_formula <- formulate(combined_inputs$out_variables, combined_inputs$anova_factors)
    formula_obj <- as.formula(user_formula)
    all_factors <- all.vars(formula_obj)[-1]
    
    # Debug: show filtered data size
    cat("DEBUG: filteredData() nrows =", nrow(filteredData()), "\n")
    if ('dbscan_cluster' %in% names(filteredData())) {
      cat("DEBUG: dbscan_cluster in filtered data: unique =", 
          length(unique(filteredData()$dbscan_cluster)),
          ", levels =", nlevels(filteredData()$dbscan_cluster), "\n")
    }
    
    cardinality <- analyzeFactorCardinality(filteredData(), all_factors)
    
    # Match model_resolver.R logic: >10 suggests lmer, but >20 will fallback to aov (unless forced)
    max_card <- if (length(cardinality) > 0) max(cardinality, na.rm = TRUE) else 0
    high_card_factors <- names(cardinality[cardinality > 10])
    blocking_candidates <- intersect(high_card_factors, c('dbscan_cluster', 'timestamp_group'))
    
    # Check if user is forcing a model type
    force_method <- if (!is.null(input$model_selection_method) && input$model_selection_method != 'auto') {
      input$model_selection_method
    } else {
      NULL
    }
    
    # Predict model type matching resolver logic (including post-fit fallback)
    use_spline <- (!is.null(force_method) && force_method == 'spline') || 
                  ('dbscan_cluster' %in% names(cardinality) && cardinality[['dbscan_cluster']] > 10)
    
    # Simulate resolver's full logic:
    # 1. Initial model selection based on cardinality
    # 2. BUT: If lmer/spline is fitted AND cardinality >20, it will fall back to aov (unless forced)
    
    if (max_card > 20) {
      # High cardinality case
      if (!is.null(force_method) && force_method %in% c('lmer', 'spline')) {
        # User forced lmer/spline - resolver will TRY but likely timeout and fallback
        model_type <- 'aov'  # Predict fallback due to timeout (realistic expectation)
        blocking_candidates <- character(0)
        use_spline <- FALSE
        force_will_fallback <- TRUE
      } else {
        # Auto mode with >20 - resolver skips lmer entirely
        model_type <- 'aov'
        blocking_candidates <- character(0)
        use_spline <- FALSE
        force_will_fallback <- FALSE
      }
    } else if (length(blocking_candidates) > 0 || (!is.null(force_method) && force_method %in% c('lmer', 'spline'))) {
      # 11-20 levels OR forced: lmer/spline will be used
      model_type <- if (!is.null(force_method) && force_method == 'spline') 'spline' else 'lmer'
      force_will_fallback <- FALSE
    } else {
      # ≤10 levels: classic aov
      model_type <- 'aov'
      force_will_fallback <- FALSE
    }
    
    # Build actual formula that will be used
    if (length(blocking_candidates) > 0) {
      # Mixed model: build lmer formula
      fixed_factors <- setdiff(all_factors, blocking_candidates)
      
      if (length(fixed_factors) == 0) {
        fixed_part <- paste(combined_inputs$out_variables, '~ 1')
      } else if (length(fixed_factors) == 1) {
        fixed_part <- paste(combined_inputs$out_variables, '~', fixed_factors[1])
      } else {
        fixed_part <- paste(combined_inputs$out_variables, '~', paste0('(', paste(fixed_factors, collapse = ' + '), ')^2'))
      }
      
      # Add random effects
      random_parts <- c()
      for (bf in blocking_candidates) {
        random_parts <- c(random_parts, paste0('(1|', bf, ')'))
      }
      
      actual_formula <- if (length(random_parts) > 0) {
        paste(fixed_part, '+', paste(random_parts, collapse = ' + '))
      } else {
        fixed_part
      }
    } else {
      # Classical ANOVA: use user formula as-is
      actual_formula <- user_formula
    }
    
    cat("DEBUG: user_formula =", user_formula, "\n")
    cat("DEBUG: actual_formula =", actual_formula, "\n")
    cat("DEBUG: cardinality =", paste(names(cardinality), cardinality, sep=':', collapse=', '), "\n")
    cat("DEBUG: model_type =", model_type, "\n")
    
    tagList(
      if (user_formula != actual_formula) {
        tags$div(
          h4(paste("User formula:", user_formula)),
          h4(paste("Actual formula:", actual_formula), style = "color: #0066cc;")
        )
      } else {
        h4(paste("Formula:", actual_formula))
      },
      if (max_card > 20 && exists('force_will_fallback') && force_will_fallback) {
        # Forced spline/lmer but will likely timeout and fallback
        tags$div(
          style = "color: #ff8800; font-size: 14px; margin-top: 5px;",
          tags$strong("⚠️ Model: "),
          sprintf("%s (likely fallback) — High cardinality (%d levels) will cause spline fit to timeout (>5s). Will fall back to classical ANOVA.", 
                  model_type, max_card)
        )
      } else if (max_card > 20 && is.null(force_method)) {
        # High cardinality warning (auto mode only)
        tags$div(
          style = "color: #ff8800; font-size: 14px; margin-top: 5px;",
          tags$strong("⚠️ Model: "),
          sprintf("%s (classical ANOVA) — High cardinality (%d levels) will skip mixed model due to performance.", 
                  model_type, max_card)
        )
      } else if (use_spline && model_type == 'spline') {
        # Spline will be used (11-20 levels)
        tags$div(
          style = "color: #9933ff; font-size: 14px; margin-top: 5px;",
          tags$strong("🌀 Model: "),
          sprintf("Spline LMM (nonlinear growth) — Using spline basis with random effects. Medium cardinality (%s: %d levels) - optimal for time-aware modeling.", 
                  paste(blocking_candidates, collapse=', '),
                  max_card)
        )
      } else if (length(blocking_candidates) > 0 || (!is.null(force_method) && force_method == 'lmer')) {
        tags$div(
          style = "color: #0066cc; font-size: 14px; margin-top: 5px;",
          tags$strong("ℹ Model: "),
          sprintf("%s (mixed model) — High-cardinality factors (%s: %s selected levels >10) moved to random effects for speed", 
                  model_type, 
                  paste(blocking_candidates, collapse=', '),
                  paste(cardinality[blocking_candidates], collapse=', '))
        )
      } else {
        tags$div(
          style = "color: #666; font-size: 14px; margin-top: 5px;",
          tags$strong("ℹ Model: "),
          paste(model_type, "(classical ANOVA)"),
          tags$span(
            style = "margin-left: 10px; color: #999;",
            sprintf("— Cardinality: %s", paste(paste(names(cardinality), cardinality, sep='='), collapse=', '))
          )
        )
      }
    )
  })
  
  # Create main plot function (reusable for download)
  createMainPlot <- function() {
    req(filteredData())
    
    data_to_plot <- withBenchmark('Shiny:PrepPlotData', filteredData())
    
    # Add group letters if requested
    if (!is.null(input$color_by) && input$color_by == 'group_letter' && 
        !is.null(analysisResults()) && !is.null(analysisResults()$tukey) && 
        analysisResults()$tukey$success) {
      
      tukey_results <- analysisResults()$tukey$results
      data_to_plot$group_letter <- NA_character_
      
      cat("DEBUG: Adding group letters. Tukey factors:", paste(names(tukey_results), collapse=", "), "\n")
      
      # Join letters from all Tukey factors
      for (factor_name in names(tukey_results)) {
        letters_df <- tukey_results[[factor_name]]$letters
        cat("DEBUG: Factor", factor_name, "- letters_df rows:", nrow(letters_df), 
            "- in data:", factor_name %in% names(data_to_plot), "\n")
        
        if (!is.null(letters_df) && nrow(letters_df) > 0 && factor_name %in% names(data_to_plot)) {
          cat("DEBUG: letters_df for", factor_name, ":\n")
          print(head(letters_df))
          cat("DEBUG: data factor levels:", paste(head(unique(data_to_plot[[factor_name]]), 10), collapse=", "), "\n")
          cat("DEBUG: data factor class:", class(data_to_plot[[factor_name]]), "\n")
          cat("DEBUG: letters group class:", class(letters_df$group), "\n")
          
          # Ensure matching types for join
          if (is.factor(data_to_plot[[factor_name]])) {
            letters_df$group <- as.character(letters_df$group)
          }
          
          # Create temporary column for this factor's letters
          temp_col <- paste0("_", factor_name, "_letter")
          data_to_plot <- data_to_plot %>%
            left_join(letters_df, by = setNames("group", factor_name)) %>%
            mutate(!!temp_col := letter) %>%
            select(-letter)
          
          cat("DEBUG: After join, temp_col", temp_col, "has", sum(!is.na(data_to_plot[[temp_col]])), "non-NA values\n")
        }
      }
      
      # Combine all letters into group_letter
      letter_cols <- names(data_to_plot)[grepl("^_[^_]+_letter$", names(data_to_plot))]
      if (length(letter_cols) > 0) {
        data_to_plot$group_letter <- apply(data_to_plot[letter_cols], 1, function(row) {
          non_na_letters <- row[!is.na(row)]
          if (length(non_na_letters) > 0) {
            paste(non_na_letters, collapse = "")
          } else {
            NA_character_
          }
        })
        
        # Clean up temporary columns
        data_to_plot <- data_to_plot %>% select(-all_of(letter_cols))
      }
      
      cat("DEBUG: Final group_letter - NAs:", sum(is.na(data_to_plot$group_letter)), "/ Total:", nrow(data_to_plot), "\n")
    }
    
    # Determine plot type and aesthetics
    use_color <- !is.null(input$color_by) && input$color_by != 'none'
    color_var <- if (use_color) input$color_by else NULL
    
    # Faceting
    facet_ok <- !is.null(input$facet_formula) && nzchar(input$facet_formula) && input$facet_formula != '~ .'
    facet_fml <- if (facet_ok) tryCatch(as.formula(input$facet_formula), error = function(e) NULL) else NULL
    
    # Create plot based on plot type
    if (!is.null(combined_inputs$timeseries_plot) && combined_inputs$timeseries_plot) {
      # Timeseries plot
      time_var <- if ("timestamp_group" %in% names(data_to_plot)) {
        "timestamp_group"
      } else if ("time_numeric" %in% names(data_to_plot)) {
        "time_numeric"
      } else {
        "dbscan_cluster"
      }
      
      p <- data_to_plot %>%
        ggplot(aes(x = !!sym(time_var), y = !!sym(combined_inputs$out_variables))) +
        stat_summary(fun = median, geom = "line", aes(group = interaction(!!!syms(if (use_color) c(color_var) else character(0))))) +
        stat_summary(fun = median, geom = "point", size = 2) +
        getStatfarmerTheme() +
        labs(
          title = paste("Timeseries of", combined_inputs$out_variables),
          x = if (time_var == "timestamp_group") "Date/Time" else if (time_var == "time_numeric") "Time (days)" else "Time Cluster",
          y = combined_inputs$out_variables
        )
      
      # Rotate x-axis labels if using timestamps (better readability)
      if (time_var == "timestamp_group") {
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      
      if (use_color) p <- p + aes(color = !!sym(color_var)) + viridis::scale_color_viridis(discrete = TRUE, end = 0.9)
      if (facet_ok) p <- p + facet_grid(facet_fml)
    } else {
      # Boxplot
      p <- data_to_plot %>%
        ggplot(aes(x = !!sym(combined_inputs$factor_grouping), y = !!sym(combined_inputs$out_variables))) +
        geom_boxplot(outlier.alpha = if (combined_inputs$outliers_plot) 0.6 else 0) +
        getStatfarmerTheme() +
        labs(
          title = paste("Distribution of", combined_inputs$out_variables),
          x = combined_inputs$factor_grouping,
          y = combined_inputs$out_variables
        )
      if (use_color) p <- p + aes(color = !!sym(color_var)) + viridis::scale_color_viridis(discrete = TRUE, end = 0.9)
      if (facet_ok) p <- p + facet_grid(facet_fml)
    }
    
    return(withBenchmark('Shiny:RenderPlot', p))
  }
  
  # Main plot
  output$distPlot <- renderPlot({
    createMainPlot()
  })
  
  # Download plot as SVG
  output$savePlot <- downloadHandler(
    filename = function() {
      paste0("statfarmer_plot_", Sys.Date(), ".svg")
    },
    content = function(file) {
      # Recreate the plot for download
      plot_obj <- createMainPlot()
      ggsave(file, plot = plot_obj, 
             width = input$plot_width, height = input$plot_height, 
             units = "mm", device = "svg")
    }
  )
  
  # Perform statistical analysis (only recomputes when submit button pressed)
  analysisResults <- reactive({
    # Only recompute when submit button is pressed (not when filters change)
    input$submit
    
    # Isolate all dependencies so they don't trigger recomputation
    data <- isolate(filteredData())
    out_var <- isolate(combined_inputs$out_variables)
    anova_fac <- isolate(combined_inputs$anova_factors)
    tukey_fac <- isolate(combined_inputs$tukey_factors)
    force_method <- isolate({
      if (!is.null(input$model_selection_method) && input$model_selection_method != 'auto') {
        input$model_selection_method
      } else {
        NULL
      }
    })
    
    req(data)
    req(out_var)
    req(anova_fac)
    
    # Create formula
    formula_text <- formulate(out_var, anova_fac)
    
    cat("DEBUG analysisResults: COMPUTING ANOVA (only on submit button press)\n")
    anova_results <- tryCatch({
      withBenchmark('Shiny:ANOVA', performANOVA(data, formula_text, out_var, force_method = force_method))
    }, error = function(e) {
      list(success = FALSE, error = paste("ANOVA failed:", e$message))
    })
    
    # Perform Tukey if ANOVA succeeded (skip for high-cardinality)
    tukey_results <- NULL
    tukey_skip_reason <- NULL
    if (anova_results$success && !is.null(tukey_fac)) {
      # Skip Tukey for ANY model with high cardinality (emmeans will be too slow)
      max_card <- if (!is.null(anova_results$cardinality)) max(anova_results$cardinality, na.rm = TRUE) else 0
      if (max_card > 20) {
        tukey_skip_reason <- paste0("Skipped: too many factor levels (", max_card, ") for pairwise comparisons. ",
                                    "Use growth summaries or effect sizes for high-cardinality factors.")
      } else {
        tukey_results <- tryCatch({
          withBenchmark('Shiny:Tukey', performTukey(anova_results$model, tukey_fac))
        }, error = function(e) {
          list(success = FALSE, error = paste("Tukey failed:", e$message))
        })
      }
    }
    
    return(list(
      anova = anova_results,
      tukey = tukey_results,
      tukey_skip_reason = tukey_skip_reason,
      formula = formula_text
    ))
  })
  
  # Raw data summary
  output$raw <- renderText({
    req(filteredData())
    data_summary <- paste(
      "Data dimensions:", nrow(filteredData()), "rows,", ncol(filteredData()), "columns\n",
      "Variables:", paste(colnames(filteredData()), collapse = ", ")
    )
    return(data_summary)
  })
  
  # Descriptive statistics
  output$DescriptiveTable <- DT::renderDataTable({
    req(filteredData())
    req(combined_inputs$out_variables)
    req(combined_inputs$anova_factors)
    
    # Calculate descriptive statistics independently of ANOVA
    tryCatch({
      data <- filteredData()
      out_var <- combined_inputs$out_variables
      factors <- combined_inputs$anova_factors
      
      cat("DEBUG DescriptiveTable: calculating desc stats for", out_var, "by", paste(factors, collapse=", "), "\n")
      cat("DEBUG DescriptiveTable: data nrows =", nrow(data), "\n")
      cat("DEBUG DescriptiveTable: testing e1071::skewness on sample data...\n")
      test_data <- c(1, 2, 3, 4, 5)
      if (requireNamespace("e1071", quietly = TRUE)) {
        test_skew <- tryCatch(e1071::skewness(test_data, type = 2), error = function(e) paste("ERROR:", e$message))
        cat("DEBUG DescriptiveTable: e1071 test result =", test_skew, "\n")
      } else {
        cat("DEBUG DescriptiveTable: e1071 not available\n")
      }
      
      # Check if all required columns exist
      if (!out_var %in% names(data)) {
        cat("DEBUG DescriptiveTable: ERROR - out_var", out_var, "not in data\n")
        return(data.frame(Message = paste("Variable", out_var, "not found in data")))
      }
      
      missing_factors <- factors[!factors %in% names(data)]
      if (length(missing_factors) > 0) {
        cat("DEBUG DescriptiveTable: ERROR - factors", paste(missing_factors, collapse=", "), "not in data\n")
        return(data.frame(Message = paste("Factors", paste(missing_factors, collapse=", "), "not found in data")))
      }
      
      # Calculate descriptive statistics
      desc_stats <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(factors))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean = mean(.data[[out_var]], na.rm = TRUE),
          median = median(.data[[out_var]], na.rm = TRUE),
          sd = sd(.data[[out_var]], na.rm = TRUE),
          se = sd / sqrt(n),
          min = min(.data[[out_var]], na.rm = TRUE),
          max = max(.data[[out_var]], na.rm = TRUE),
          q25 = quantile(.data[[out_var]], 0.25, na.rm = TRUE),
          q75 = quantile(.data[[out_var]], 0.75, na.rm = TRUE),
          iqr = q75 - q25,
          cv = (sd / mean) * 100,  # Coefficient of variation (%)
          .groups = 'drop'
        )
      
      # Add skewness and kurtosis using a simpler approach
      desc_stats$skewness <- NA_real_
      desc_stats$kurtosis <- NA_real_
      
      # Calculate skewness and kurtosis for each group
      for (i in 1:nrow(desc_stats)) {
        # Get the group values
        group_data <- data
        for (factor in factors) {
          group_data <- group_data[group_data[[factor]] == desc_stats[[factor]][i], ]
        }
        
        x <- group_data[[out_var]][!is.na(group_data[[out_var]])]
        
        if (length(x) > 2) {
          desc_stats$skewness[i] <- e1071::skewness(x, type = 2)
        }
        if (length(x) > 3) {
          desc_stats$kurtosis[i] <- e1071::kurtosis(x, type = 2)
        }
      }
      
      cat("DEBUG DescriptiveTable: calculated", nrow(desc_stats), "rows of desc stats\n")
      cat("DEBUG DescriptiveTable: checking packages - e1071:", requireNamespace("e1071", quietly = TRUE), 
          ", moments:", requireNamespace("moments", quietly = TRUE), "\n")
      
      # Apply precision formatting
      desc_stats <- desc_stats %>%
        mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 3)))
      
      return(desc_stats)
    }, error = function(e) {
      cat("DEBUG DescriptiveTable: ERROR =", e$message, "\n")
      return(data.frame(Message = paste("Error calculating descriptive statistics:", e$message)))
    })
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # ANOVA results
  output$anovaTable <- DT::renderDataTable({
    req(analysisResults())
    if (analysisResults()$anova$success) {
      anova_table <- analysisResults()$anova$anova_table %>%
        mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 3)))
      return(anova_table)
    } else {
      return(data.frame(Message = paste("ANOVA Error:", analysisResults()$anova$error)))
    }
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # ANOVA model info
  output$anovaResults <- renderText({
    cat("DEBUG anovaResults: renderText called\n")
    req(analysisResults())
    cat("DEBUG anovaResults: analysisResults() is available\n")
    if (analysisResults()$anova$success) {
      model_type <- if (!is.null(analysisResults()$anova$model_type)) {
        analysisResults()$anova$model_type
      } else {
        'aov'
      }
      
      info_lines <- c(
        paste("Model type:", model_type),
        paste("Formula:", analysisResults()$formula)
      )
      
      cat("DEBUG anovaResults: analysisResults()$anova$cardinality =", 
          if (is.null(analysisResults()$anova$cardinality)) "NULL" else paste(names(analysisResults()$anova$cardinality), analysisResults()$anova$cardinality, sep=':', collapse=', '), "\n")
      cat("DEBUG anovaResults: analysisResults() structure =", paste(names(analysisResults()), collapse=", "), "\n")
      cat("DEBUG anovaResults: analysisResults()$anova structure =", paste(names(analysisResults()$anova), collapse=", "), "\n")
      if (!is.null(analysisResults()$anova$cardinality)) {
        card_str <- paste(names(analysisResults()$anova$cardinality), 
                         analysisResults()$anova$cardinality, 
                         sep=':', collapse=', ')
        info_lines <- c(info_lines, paste("Factor levels:", card_str))
      }
      
      if (!is.null(analysisResults()$anova$blocking_factors) && 
          length(analysisResults()$anova$blocking_factors) > 0) {
        info_lines <- c(info_lines, 
                       paste("Random effects:", paste(analysisResults()$anova$blocking_factors, collapse=', ')))
      }
      
      paste(info_lines, collapse = '\n')
    } else {
      paste("ANOVA failed:", analysisResults()$anova$error)
    }
  })
  
  # Tukey results
  output$tukeyTable <- DT::renderDataTable({
    req(analysisResults())
    if (!is.null(analysisResults()$tukey_skip_reason)) {
      return(data.frame(Message = analysisResults()$tukey_skip_reason))
    } else if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
      tukey_results <- analysisResults()$tukey$results
      # Combine all Tukey results into one table
      all_results <- do.call(rbind, lapply(tukey_results, function(x) x$results)) %>%
        mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 3)))
      return(all_results)
    } else {
      return(data.frame(Message = "No Tukey results available"))
    }
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Group letters
  output$tukeyLetters <- DT::renderDataTable({
    req(analysisResults())
    if (!is.null(analysisResults()$tukey_skip_reason)) {
      return(data.frame(Message = analysisResults()$tukey_skip_reason))
    } else if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
      tukey_results <- analysisResults()$tukey$results
      cat("DEBUG tukeyLetters: number of factors =", length(tukey_results), "\n")
      cat("DEBUG tukeyLetters: factor names =", paste(names(tukey_results), collapse=", "), "\n")
      
      # Combine all group letters into one table
      all_letters <- data.frame()
      for (factor_name in names(tukey_results)) {
        factor_letters <- tukey_results[[factor_name]]$letters
        cat("DEBUG tukeyLetters: factor", factor_name, "letters is.null =", is.null(factor_letters), 
            "nrow =", if (!is.null(factor_letters)) nrow(factor_letters) else NA, "\n")
        
        if (!is.null(factor_letters)) {
          cat("DEBUG tukeyLetters: factor", factor_name, "letters structure:\n")
          print(str(factor_letters))
        }
        
        if (!is.null(factor_letters) && nrow(factor_letters) > 0) {
          factor_letters$factor <- factor_name
          all_letters <- rbind(all_letters, factor_letters)
        }
      }
      
      cat("DEBUG tukeyLetters: total rows in all_letters =", nrow(all_letters), "\n")
      
      if (nrow(all_letters) > 0) {
        return(all_letters)
      } else {
        return(data.frame(Message = "No group letters available"))
      }
    } else {
      return(data.frame(Message = "No group letters available"))
    }
  }, options = list(pageLength = 10, scrollX = TRUE))

  # Assumption checks
  output$assumptionsTable <- DT::renderDataTable({
    req(analysisResults())
    if (analysisResults()$anova$success) {
      # derive first factor on RHS for variance tests
      rhs_terms <- all.vars(as.formula(analysisResults()$formula))[-1]
      main_factor <- if (length(rhs_terms) > 0) rhs_terms[1] else NULL
      checks <- performAssumptionChecks(
        analysisResults()$anova$model,
        filteredData(),
        combined_inputs$out_variables,
        main_factor
      )
      checks <- checks %>% mutate(p_value = signif(p_value, 3))
      return(checks)
    } else {
      return(data.frame(Message = "No ANOVA model available"))
    }
  }, options = list(dom = 't', paging = FALSE))
  
  # ANOVA interpretation
  output$anovaResults <- renderText({
    req(analysisResults())
    if (analysisResults()$anova$success) {
      anova_table <- analysisResults()$anova$anova_table
      significant_factors <- anova_table[anova_table$p.value < 0.05 & !is.na(anova_table$p.value), ]
      
      if (nrow(significant_factors) > 0) {
        interpretation <- paste(
          "Significant factors (p < 0.05):\n",
          paste(significant_factors$term, collapse = ", "), "\n\n",
          "Model R-squared:", round(1 - (anova_table$sumsq[nrow(anova_table)] / sum(anova_table$sumsq)), 3)
        )
      } else {
        interpretation <- "No significant factors found (p < 0.05)"
      }
      return(interpretation)
    } else {
      return("Cannot interpret ANOVA results due to errors")
    }
  })
  
  # Placeholder UI outputs
  output$RAW <- renderUI(h4("Raw Data"))
  output$DESCvarNames <- renderUI(h4("Descriptive Statistics"))
  output$ANOVAvarNames <- renderUI(h4("ANOVA Results"))
  output$TUKEYvarNames <- renderUI(h4("Tukey HSD Results"))
  output$TUKEYLvarNames <- renderUI(h4("Group Letters"))
  
  # Model info box
  output$model_info_box <- renderUI({
    req(analysisResults())
    if (analysisResults()$anova$success) {
      model_type <- analysisResults()$anova$model_type
      cardinality <- analysisResults()$anova$cardinality
      blocking_factors <- analysisResults()$anova$blocking_factors
      model <- analysisResults()$anova$model
      
      # Build model description
      type_desc <- switch(model_type,
        "aov" = "Classic ANOVA (fixed effects)",
        "lmer" = "Linear Mixed Model (random effects)",
        "lm" = "Linear Model (fixed effects, spline fallback)",
        model_type
      )
      
      # Get actual formula from fitted model
      actual_formula <- tryCatch({
        if (inherits(model, 'lmerMod') || inherits(model, 'merMod')) {
          deparse(formula(model))
        } else {
          deparse(formula(model))
        }
      }, error = function(e) {
        analysisResults()$formula
      })
      
      # Build info text
      info_lines <- c(
        paste0("<strong>Model Type:</strong> ", type_desc),
        paste0("<strong>Actual Formula:</strong> <code>", actual_formula, "</code>"),
        paste0("<strong>Factor Cardinality:</strong> ", paste(names(cardinality), "=", cardinality, collapse=", "))
      )
      
      # Warn about high cardinality fallback
      max_card <- if (length(cardinality) > 0) max(cardinality, na.rm = TRUE) else 0
      if (max_card > 20 && model_type == 'aov') {
        info_lines <- c(info_lines, 
          paste0("<strong>⚠️ Note:</strong> High cardinality (", max_card, " levels) detected. ",
                 "Fell back to ANOVA for performance. Consider using growth summaries for time analysis."))
      }
      
      if (length(blocking_factors) > 0) {
        info_lines <- c(info_lines, paste0("<strong>Blocking Factors:</strong> ", paste(blocking_factors, collapse=", ")))
      }
      
      # Check if spline was used (look for time_numeric in model)
      model_frame_names <- tryCatch({
        if (inherits(model, 'lmerMod') || inherits(model, 'merMod')) {
          names(model@frame)
        } else {
          names(model$model)
        }
      }, error = function(e) character(0))
      
      if (model_type %in% c("lmer", "lm") && "time_numeric" %in% model_frame_names) {
        # Extract spline details from model
        spline_info <- "<strong>Spline Basis:</strong> Used (nonlinear growth modeling)"
        
        # Try to detect number of knots from formula/model terms
        formula_str <- tryCatch(deparse(formula(model)), error = function(e) "")
        k_match <- regexpr("k\\s*=\\s*([0-9]+)", formula_str)
        if (k_match > 0) {
          k_val <- regmatches(formula_str, k_match)
          k_num <- gsub("k\\s*=\\s*", "", k_val)
          spline_info <- paste0(spline_info, " with ", k_num, " knots (adaptive)")
        }
        
        # Check for random slopes in spline
        if (inherits(model, 'lmerMod') && grepl("bs\\(time_numeric.*\\|", formula_str)) {
          spline_info <- paste0(spline_info, " + random slopes")
        }
        
        info_lines <- c(info_lines, spline_info)
      }
      
      # Show random effects structure if lmer/lm with random effects
      if (model_type == "lmer" && inherits(model, 'lmerMod')) {
        re_names <- names(lme4::ranef(model))
        if (length(re_names) > 0) {
          info_lines <- c(info_lines, 
            paste0("<strong>Random Effects:</strong> ", paste(paste0("(1|", re_names, ")"), collapse=" + ")))
        }
      }
      
      # Color code based on model type
      box_class <- switch(model_type,
        "aov" = if (max_card > 20) "warning" else "info",
        "lmer" = "primary",
        "lm" = "warning",
        "info"
      )
      
      div(
        class = paste0("alert alert-", box_class),
        style = "margin-top: 10px;",
        HTML(paste(info_lines, collapse = "<br/>"))
      )
    }
  })
  
  # Raw data table
  output$raw_flex <- DT::renderDT({
    req(filteredData())
    DT::datatable(
      filteredData(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = "400px"
      )
    )
  })
  
  # Descriptive statistics table
  output$descriptiveTable_flex <- DT::renderDT({
    req(filteredData())
    req(combined_inputs$out_variables)
    req(combined_inputs$anova_factors)
    
    # Calculate descriptive statistics independently of ANOVA
    tryCatch({
      data <- filteredData()
      out_var <- combined_inputs$out_variables
      factors <- combined_inputs$anova_factors
      
      # Check if all required columns exist
      if (!out_var %in% names(data)) {
        return(DT::datatable(data.frame(Message = paste("Variable", out_var, "not found in data"))))
      }
      
      missing_factors <- factors[!factors %in% names(data)]
      if (length(missing_factors) > 0) {
        return(DT::datatable(data.frame(Message = paste("Factors", paste(missing_factors, collapse=", "), "not found in data"))))
      }
      
      # Calculate descriptive statistics
      desc_stats <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(factors))) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean = mean(.data[[out_var]], na.rm = TRUE),
          median = median(.data[[out_var]], na.rm = TRUE),
          sd = sd(.data[[out_var]], na.rm = TRUE),
          se = sd / sqrt(n),
          min = min(.data[[out_var]], na.rm = TRUE),
          max = max(.data[[out_var]], na.rm = TRUE),
          q25 = quantile(.data[[out_var]], 0.25, na.rm = TRUE),
          q75 = quantile(.data[[out_var]], 0.75, na.rm = TRUE),
          iqr = q75 - q25,
          cv = (sd / mean) * 100,  # Coefficient of variation (%)
          .groups = 'drop'
        )
      
      # Add skewness and kurtosis using a simpler approach
      desc_stats$skewness <- NA_real_
      desc_stats$kurtosis <- NA_real_
      
      # Calculate skewness and kurtosis for each group
      for (i in 1:nrow(desc_stats)) {
        # Get the group values
        group_data <- data
        for (factor in factors) {
          group_data <- group_data[group_data[[factor]] == desc_stats[[factor]][i], ]
        }
        
        x <- group_data[[out_var]][!is.na(group_data[[out_var]])]
        
        if (length(x) > 2) {
          desc_stats$skewness[i] <- e1071::skewness(x, type = 2)
        }
        if (length(x) > 3) {
          desc_stats$kurtosis[i] <- e1071::kurtosis(x, type = 2)
        }
      }
      
      # Apply precision formatting
      desc_stats <- desc_stats %>%
        mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(., 3)))
      
      DT::datatable(
        desc_stats,
        options = list(
          pageLength = 15,
          scrollX = TRUE
        )
      )
    }, error = function(e) {
      DT::datatable(data.frame(Message = paste("Error calculating descriptive statistics:", e$message)))
    })
  })
  
  # ANOVA table
  output$anovaTable_flex <- DT::renderDT({
    req(analysisResults())
    if (analysisResults()$anova$success) {
      DT::datatable(
        analysisResults()$anova$anova_table,
        options = list(
          pageLength = 15,
          scrollX = TRUE
        )
      )
    } else {
      DT::datatable(data.frame(Message = "Error calculating ANOVA"))
    }
  })
  
  # Tukey table
  output$tukeyTable_flex <- DT::renderDT({
    req(analysisResults())
    if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
      # Combine all Tukey results
      all_results <- data.frame()
      for (factor_name in names(analysisResults()$tukey$results)) {
        factor_results <- analysisResults()$tukey$results[[factor_name]]$results
        factor_results$factor <- factor_name
        all_results <- rbind(all_results, factor_results)
      }
      
      if (nrow(all_results) > 0) {
        DT::datatable(
          all_results,
          options = list(
            pageLength = 15,
            scrollX = TRUE
          )
        )
      } else {
        DT::datatable(data.frame(Message = "No Tukey results available"))
      }
    } else {
      DT::datatable(data.frame(Message = "No Tukey results available"))
    }
  })
  
  # Group letters table
  output$tukeyLetters_flex <- DT::renderDT({
    req(analysisResults())
    if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
      # Combine all group letters
      all_letters <- data.frame()
      for (factor_name in names(analysisResults()$tukey$results)) {
        factor_letters <- analysisResults()$tukey$results[[factor_name]]$letters
        factor_letters$factor <- factor_name
        all_letters <- rbind(all_letters, factor_letters)
      }
      
      if (nrow(all_letters) > 0) {
        DT::datatable(
          all_letters,
          options = list(
            pageLength = 15,
            scrollX = TRUE
          )
        )
      } else {
        DT::datatable(data.frame(Message = "No group letters available"))
      }
    } else {
      DT::datatable(data.frame(Message = "No group letters available"))
    }
  })
  
  # ANOVA plot
  output$ANOVAPlot <- renderPlot({
    req(analysisResults())
    req(combined_inputs$factor_grouping)
    
    if (analysisResults()$anova$success) {
      createANOVAPlot(
        filteredData(),
        analysisResults()$formula,
        combined_inputs$out_variables,
        combined_inputs$factor_grouping
      )
    } else {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("Error creating ANOVA plot:", analysisResults()$anova$error), 
                size = 6) +
        getStatfarmerTheme() +
        theme_void()
    }
  })
  
  return(list(
    analysisResults = analysisResults
  ))
}
