# StatFaRmer Reactive UI Management
# UI updates and input management

#' Create reactive UI management functions
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param combined_inputs Reactive values for combined inputs
#' @param projectUISelections Reactive function for UI selections
#' @param projectData Reactive function for project data
#' @return List of reactive functions
createReactiveUI <- function(input, output, session, combined_inputs, projectUISelections, projectData) {
  
  # Update UI selections when project changes
  updateProjectSelections <- function() {
    observeEvent(input$selected_project, {
      selections <- withBenchmark('Shiny:UI:UpdateSelections', projectUISelections())
      logEvent('INFO', 'ui.project_changed', list(project = input$selected_project))
      
      # Update ANOVA factors (preserve current user order if available)
      anova_choices <- setdiff(selections$anova_factors, 'unit')
      current_anova <- isolate(input$anova_factors)
      anova_selected <- intersect(current_anova, anova_choices)
      if (length(anova_selected) == 0) {
        # Fallback to a sensible default but keep displayed order from choices
        anova_selected <- intersect(anova_choices, c('treatment','dbscan_cluster'))
      }
      updateSelectInput(session, 'anova_factors', choices = anova_choices, selected = anova_selected)

      # Update Tukey factors (preserve current user order if available)
      tukey_choices <- setdiff(selections$anova_factors, 'unit')
      current_tukey <- isolate(input$tukey_factors)
      tukey_selected <- intersect(current_tukey, tukey_choices)
      if (length(tukey_selected) == 0) {
        tukey_selected <- intersect(tukey_choices, c('treatment','dbscan_cluster'))
      }
      updateSelectInput(session, 'tukey_factors', choices = tukey_choices, selected = tukey_selected)
      
      # Update color by (hide Group Letter)
      updateSelectInput(session, 'color_by', 
                       choices = setdiff(selections$anova_factors, 'unit'), #'Group Letter' = 'group_letter' is hidden
                       selected = if ('treatment' %in% selections$anova_factors) 'treatment' else setdiff(selections$anova_factors, 'unit')[1])
      
      # Update grouping factor
      updateSelectInput(session, 'factor_grouping', 
                       choices = setdiff(selections$anova_factors, 'unit'),
                       selected = if ('treatment' %in% selections$anova_factors) 'treatment' else setdiff(selections$anova_factors, 'unit')[1])
      
      # Update treatments
      updateSelectInput(session, 'treatments', 
                       choices = selections$treatments,
                       selected = selections$treatments)
      
      # Update cultivars (default: all)
      updateSelectInput(session, 'cultivars', 
                       choices = selections$cultivars,
                       selected = selections$cultivars)
      
      # Update timestamp groups (first, middle, last)
      updatePickerInput(session, 'timestamp_groups', 
                       choices = selections$named_timestamp_groups,
                       selected = {
                         n <- length(selections$timestamp_groups)
                         if (n >= 3) {
                           as.character(selections$timestamp_groups[c(1, ceiling(n/2), n)])
                         } else {
                           as.character(selections$timestamp_groups)
                         }
                       })
      
      # Update output variables
      updateSelectInput(session, 'out_variables', 
                       choices = selections$out_variables,
                       selected = selections$out_variables[1])
      
      # Reset combined inputs
      combined_inputs$factor_grouping <- NULL
      combined_inputs$factor_levels <- NULL
      combined_inputs$treatments <- NULL
      combined_inputs$cultivars <- NULL
      combined_inputs$timestamp_groups <- NULL
      combined_inputs$out_variables <- NULL
      combined_inputs$anova_factors <- NULL
      combined_inputs$tukey_factors <- NULL
    }, priority = 100)
  }
  
  # Update factor levels when grouping factor or project changes
  updateFactorLevels <- function(projectData) {
    observeEvent(list(input$factor_grouping, input$selected_project), {
      current_data <- projectData()
      if (nrow(current_data$merged_table) > 0) {
        unique_groups <- getUnique(current_data$merged_table, input$factor_grouping)
        
        if (!setequal(unique_groups, input$factor_levels)) {
          isolate(
            updateSelectizeInput(
              session,
              'factor_levels',
              choices = c(),
              selected = c(),
              server = TRUE
            )
          )
          updateSelectizeInput(
            session,
            'factor_levels',
            choices = unique_groups,
            selected = unique_groups,
            server = TRUE
          )
        }
      }
    }, ignoreInit = FALSE, priority = 200)
  }
  
  # Update cultivars when factor levels change
  updateCultivars <- function(projectData) {
    observeEvent(c(input$factor_levels, input$factor_grouping), {
      # Avoid clearing cultivars while factor_levels is being reset during project change
      req(!is.null(input$factor_grouping))
      req(length(input$factor_levels) > 0)
      current_data <- projectData()
      if (nrow(current_data$merged_table) > 0) {
        unique_cultivars <- current_data$merged_table %>%
          filter(!!sym(input$factor_grouping) %in% input$factor_levels) %>%
          dplyr::select(dplyr::all_of('cultivar')) %>%
          dplyr::distinct() %>%
          dplyr::arrange(cultivar) %>%
          dplyr::pull(1)
        
        if (!setequal(unique_cultivars, input$cultivars)) {
          isolate(
            updateSelectizeInput(
              session,
              'cultivars',
              choices = c(),
              selected = c(),
              server = TRUE
            )
          )
          updateSelectizeInput(
            session,
            'cultivars',
            choices = unique_cultivars,
            selected = unique_cultivars,
            server = TRUE
          )
        }
      }
    }, ignoreInit = FALSE, priority = 100)
  }
  
  # Create initial combined inputs
  createInitialCombinedInputs <- function() {
    reactive({
      list(
        factor_grouping = isolate(input$factor_grouping),
        factor_levels = isolate(input$factor_levels),
        treatments = isolate(input$treatments),
        cultivars = isolate(input$cultivars),
        timestamp_groups = {
          # Convert character timestamp_groups back to POSIXct values from the actual data
          selected_timestamp_chars <- isolate(unname(input$timestamp_groups))
          if (length(selected_timestamp_chars) > 0) {
            # Get the actual POSIXct values from the data that match the selected character strings
            current_data <- projectData()
            if (nrow(current_data$merged_table) > 0) {
              actual_timestamp_groups <- sort(unique(current_data$merged_table$timestamp_group))
              # Find matching POSIXct values for the selected character strings
              actual_timestamp_groups[
                as.character(actual_timestamp_groups) %in% selected_timestamp_chars
              ]
            } else {
              as.POSIXct(selected_timestamp_chars, tz = "UTC")
            }
          } else {
            as.POSIXct(character(0))
          }
        },
        out_variables = isolate(input$out_variables),
        facet_formula = isolate(as.formula(input$facet_formula)),
        anova_factors = isolate(input$anova_factors),
        tukey_factors = isolate(input$tukey_factors),
        timeseries_plot = isolate(input$timeseries_plot),
        outliers_plot = isolate(input$outliers_plot)
      )
    })
  }
  
  # Update combined inputs on initialization
  updateCombinedInputsInit <- function() {
    # Use observeEvent with once=TRUE to wait for first non-NULL value
    observeEvent(input$timestamp_groups, {
      cat("DEBUG updateCombinedInputsInit: timestamp_groups length =", length(input$timestamp_groups), "\n")
      
      initial_combined_inputs <- createInitialCombinedInputs()()
      
      combined_inputs$factor_grouping <- initial_combined_inputs$factor_grouping
      combined_inputs$factor_levels <- initial_combined_inputs$factor_levels
      combined_inputs$treatments <- initial_combined_inputs$treatments
      combined_inputs$cultivars <- initial_combined_inputs$cultivars
      combined_inputs$timestamp_groups <- initial_combined_inputs$timestamp_groups
      combined_inputs$out_variables <- initial_combined_inputs$out_variables
      combined_inputs$facet_formula <- initial_combined_inputs$facet_formula
      combined_inputs$anova_factors <- initial_combined_inputs$anova_factors
      combined_inputs$tukey_factors <- initial_combined_inputs$tukey_factors
      combined_inputs$timeseries_plot <- initial_combined_inputs$timeseries_plot
      combined_inputs$outliers_plot <- initial_combined_inputs$outliers_plot
      
      cat("DEBUG updateCombinedInputsInit: combined_inputs$timestamp_groups length =", 
          length(combined_inputs$timestamp_groups), "\n")
    }, once = TRUE, ignoreNULL = TRUE, priority = 50)
  }
  
  # Update combined inputs on submit
  updateCombinedInputsSubmit <- function() {
    observeEvent(input$submit, {
      logEvent('INFO', 'ui.submit', list(
        factor_grouping = input$factor_grouping,
        factor_levels = input$factor_levels,
        treatments = input$treatments,
        cultivars = input$cultivars,
        timestamp_groups = input$timestamp_groups,
        out_variables = input$out_variables,
        facet_formula = input$facet_formula,
        anova_factors = input$anova_factors,
        tukey_factors = input$tukey_factors,
        timeseries_plot = input$timeseries_plot,
        outliers_plot = input$outliers_plot
      ))
      combined_inputs$factor_grouping = isolate(input$factor_grouping)
      combined_inputs$factor_levels = isolate(input$factor_levels)
      combined_inputs$treatments = isolate(input$treatments)
      combined_inputs$cultivars = isolate(input$cultivars)
      # Convert character timestamp_groups back to POSIXct values from the actual data
      # to prevent precision mismatch in filtering
      selected_timestamp_chars <- isolate(unname(input$timestamp_groups))
      if (length(selected_timestamp_chars) > 0) {
        # Get the actual POSIXct values from the data that match the selected character strings
        current_data <- projectData()
        if (nrow(current_data$merged_table) > 0) {
          actual_timestamp_groups <- sort(unique(current_data$merged_table$timestamp_group))
          # Find matching POSIXct values for the selected character strings
          selected_timestamp_groups <- actual_timestamp_groups[
            as.character(actual_timestamp_groups) %in% selected_timestamp_chars
          ]
          combined_inputs$timestamp_groups <- selected_timestamp_groups
        } else {
          combined_inputs$timestamp_groups <- as.POSIXct(selected_timestamp_chars, tz = "UTC")
        }
      } else {
        combined_inputs$timestamp_groups <- as.POSIXct(character(0))
      }
      combined_inputs$out_variables = isolate(input$out_variables)
      combined_inputs$facet_formula = isolate(as.formula(input$facet_formula))
      combined_inputs$anova_factors = isolate(input$anova_factors)
      combined_inputs$tukey_factors = isolate(input$tukey_factors)
      combined_inputs$timeseries_plot = isolate(input$timeseries_plot)
      combined_inputs$outliers_plot = isolate(input$outliers_plot)
    }, priority = 50)
  }
  
  return(list(
    updateProjectSelections = updateProjectSelections,
    updateFactorLevels = updateFactorLevels,
    updateCultivars = updateCultivars,
    createInitialCombinedInputs = createInitialCombinedInputs,
    updateCombinedInputsInit = updateCombinedInputsInit,
    updateCombinedInputsSubmit = updateCombinedInputsSubmit
  ))
}
