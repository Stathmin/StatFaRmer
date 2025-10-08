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
        # Return empty data structure
        list(
          merged_table = data.frame(),
          vector_of_groups = character(0)
        )
      })
    })
  })
  
  # Reactive UI selections based on current project data
  projectUISelections <- reactive({
    data <- projectData()
    if (nrow(data$merged_table) > 0) {
      prepareUISelections(data$merged_table)
    } else {
      list(
        anova_factors = character(0),
        treatments = character(0),
        cultivars = character(0),
        timestamp_groups = character(0),
        named_timestamp_groups = character(0),
        out_variables = character(0)
      )
    }
  })
  
  # Filter data based on inputs
  filteredData <- reactive({
    req(combined_inputs$factor_grouping)
    req(combined_inputs$factor_levels)
    req(combined_inputs$treatments)
    req(combined_inputs$cultivars)
    req(combined_inputs$timestamp_groups)
    req(combined_inputs$out_variables)
    
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

    filtered <- withBenchmark('Shiny:FilterData', {
      current_data$merged_table %>%
      {
        # Apply timestamp filter first if we have mapped clusters
        tmp <- if (length(mapped_clusters) > 0) {
          dplyr::filter(., dbscan_cluster %in% mapped_clusters)
        } else {
          .
        }
        
        # Then apply other filters
        tmp <- dplyr::filter(
          tmp,
          !!sym(combined_inputs$factor_grouping) %in% combined_inputs$factor_levels,
          treatment %in% combined_inputs$treatments,
          cultivar %in% combined_inputs$cultivars
        )
        
        # Drop unused factor levels early to keep downstream cardinality correct
        tmp <- droplevels(tmp)
        tmp
      }
    })
    
    # Handle outliers
    if (!combined_inputs$outliers_plot) {
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
