# StatFaRmer Reactive Data Management
# Project data loading and UI selections

#' Create reactive data loading functions
#' @param input Shiny input object
#' @param session Shiny session object
#' @param combined_inputs Reactive values list from server
#' @return List of reactive functions
createReactiveData <- function(input, session, combined_inputs) {
  
  # Reactive data loading based on selected project
  projectData <- reactive({
    req(input$selected_project)
    withBenchmark(sprintf('Shiny:LoadProject:%s', input$selected_project), {
      tryCatch({
        res <- loadProjectData(input$selected_project)
        logInfo(paste('Loaded project', input$selected_project))
        res
      }, error = function(e) {
        logError(paste('Load project failed', input$selected_project, e$message))
        # Show error message to user
        showNotification(
          paste("Error loading project", input$selected_project, ":", e$message),
          type = "error",
          duration = 10
        )
        emptyProjectData()
      })
    })
  })
  
  # Reactive UI selections based on current project data
  projectUISelections <- reactive({
    data <- projectData()
    if (nrow(data$merged_table) > 0) {
      prepareUISelections(data$merged_table)
    } else {
      emptyUISelections()
    }
  })
  
  # Filter data based on inputs
  filteredData <- reactive({
    # Trait + times are required; grouping/treatment may be soft (#8 single treatment, #10 none)
    req(combined_inputs$out_variables)
    req(combined_inputs$timestamp_groups)
    
    # Get current project data
    current_data <- projectData()
    if (nrow(current_data$merged_table) == 0) {
      return(data.frame())
    }
    
    # Filter data
    # Map selected timestamp_group -> dbscan_cluster; if mapping empty, skip this filter
    cat("DEBUG filteredData: timestamp_groups selected =", length(combined_inputs$timestamp_groups), "\n")
    cat("DEBUG filteredData: timestamp_groups class =", class(combined_inputs$timestamp_groups), "\n")
    if (length(combined_inputs$timestamp_groups) > 0) {
      cat("DEBUG filteredData: first 3 timestamp_groups =", as.character(head(combined_inputs$timestamp_groups, 3)), "\n")
    }
    
    # Convert timestamp_groups to POSIXct if needed (they might be character from UI)
    timestamp_selection <- if (is.character(combined_inputs$timestamp_groups)) {
      as.POSIXct(combined_inputs$timestamp_groups, tz = "UTC")
    } else {
      combined_inputs$timestamp_groups
    }
    
    mapped_clusters <- tryCatch({
      filtered_by_time <- current_data$merged_table %>% 
        dplyr::filter(timestamp_group %in% timestamp_selection)
      cat("DEBUG filteredData: rows after timestamp filter =", nrow(filtered_by_time), "\n")
      clusters <- getUnique(filtered_by_time, 'dbscan_cluster')
      cat("DEBUG filteredData: mapped_clusters =", length(clusters), "clusters:", paste(head(clusters, 10), collapse=','), "\n")
      clusters
    }, error = function(e) {
      cat("DEBUG filteredData: ERROR in mapping:", e$message, "\n")
      character(0)
    })

    grouping <- combined_inputs$factor_grouping
    use_grouping <- !is.null(grouping) && nzchar(grouping) && !identical(grouping, "none") &&
      grouping %in% names(current_data$merged_table)
    factor_levels <- combined_inputs$factor_levels
    treatments <- combined_inputs$treatments
    cultivars <- combined_inputs$cultivars

    filtered <- withBenchmark('Shiny:FilterData', {
      current_data$merged_table %>%
      {
        # Apply timestamp filter first if we have mapped clusters
        tmp <- if (length(mapped_clusters) > 0) {
          dplyr::filter(., dbscan_cluster %in% mapped_clusters)
        } else {
          .
        }
        
        # Grouping levels — skip when "none" or levels empty (#8/#10)
        if (use_grouping && length(factor_levels) > 0) {
          tmp <- dplyr::filter(tmp, !!sym(grouping) %in% factor_levels)
        }
        # Treatment / cultivar — skip filter when selector empty so 1-level projects still plot
        if (length(treatments) > 0 && "treatment" %in% names(tmp)) {
          tmp <- dplyr::filter(tmp, treatment %in% treatments)
        }
        if (length(cultivars) > 0 && "cultivar" %in% names(tmp)) {
          tmp <- dplyr::filter(tmp, cultivar %in% cultivars)
        }
        
        # Drop unused factor levels early to keep downstream cardinality correct
        tmp <- droplevels(tmp)
        tmp
      }
    })
    
    # Handle outliers
    if (!isTRUE(combined_inputs$outliers_plot) && "outlier" %in% names(filtered)) {
      filtered <- filtered %>% filter(!outlier)
    }
    
    return(filtered)
  })
  
  return(list(
    projectData = projectData,
    projectUISelections = projectUISelections,
    filteredData = filteredData
  ))
}
