# StatFaRmer Reactive Download Management
# Download handlers

#' Create reactive download functions
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param combined_inputs Reactive values for combined inputs
#' @param filteredData Reactive function for filtered data
#' @param analysisResults Reactive function for analysis results
#' @return List of reactive functions
createReactiveDownloads <- function(input, output, session, combined_inputs, filteredData, analysisResults) {
  
  # Download raw data (matching UI button ID)
  output$raw_flex_downloadData <- downloadHandler(
    filename = function() { 
      paste0("raw_data_", combined_inputs$out_variables, "_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      write.csv(filteredData(), file, row.names = FALSE) 
    }
  )
  
  # Download descriptive statistics (matching UI button ID)
  output$descriptiveTable_downloadData <- downloadHandler(
    filename = function() { 
      paste0("descriptive_stats_", combined_inputs$out_variables, "_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      if (analysisResults()$anova$success) {
        write.csv(analysisResults()$anova$desc_stats, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "Error calculating descriptive statistics"), file)
      }
    }
  )
  
  # Download ANOVA results (matching UI button ID)
  output$anovaTable_downloadData <- downloadHandler(
    filename = function() { 
      paste0("anova_results_", combined_inputs$out_variables, "_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      if (analysisResults()$anova$success) {
        write.csv(analysisResults()$anova$anova_table, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "Error calculating ANOVA"), file)
      }
    }
  )
  
  # Download Tukey results (matching UI button ID)
  output$tukeyTable_downloadData <- downloadHandler(
    filename = function() { 
      paste0("tukey_results_", combined_inputs$out_variables, "_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
        # Combine all Tukey results
        all_results <- data.frame()
        for (factor_name in names(analysisResults()$tukey$results)) {
          factor_results <- analysisResults()$tukey$results[[factor_name]]$results
          factor_results$factor <- factor_name
          all_results <- rbind(all_results, factor_results)
        }
        write.csv(all_results, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No Tukey results available"), file)
      }
    }
  )
  
  # Download group letters (matching UI button ID)
  output$tukeyLetters_downloadData <- downloadHandler(
    filename = function() { 
      paste0("group_letters_", combined_inputs$out_variables, "_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      if (!is.null(analysisResults()$tukey) && analysisResults()$tukey$success) {
        # Combine all group letters
        all_letters <- data.frame()
        for (factor_name in names(analysisResults()$tukey$results)) {
          factor_letters <- analysisResults()$tukey$results[[factor_name]]$letters
          factor_letters$factor <- factor_name
          all_letters <- rbind(all_letters, factor_letters)
        }
        write.csv(all_letters, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No group letters available"), file)
      }
    }
  )
  
  return(list())
}
