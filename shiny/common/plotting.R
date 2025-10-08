# Plotting utilities

library(ggplot2)
library(rlang)
library(cowplot)

#' Apply consistent cowplot theme with enhanced faceting
#' @return Theme object with cowplot base and enhanced styling
getStatfarmerTheme <- function() {
  cowplot::theme_cowplot() +
    theme(
      # Enhanced faceting control
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
      strip.text.x = element_text(margin = margin(t = 8, b = 8)),
      strip.text.y = element_text(margin = margin(l = 8, r = 8)),
      
      # Typography
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      
      # Legends
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      
      # Grid and background
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
}

#' Create ANOVA plot
createANOVAPlot <- function(data, formula, out_variable, factor_grouping) {
  tryCatch({
    p <- data %>%
      ggplot(aes(x = !!sym(factor_grouping), y = !!sym(out_variable))) +
      geom_boxplot() +
      getStatfarmerTheme() +
      labs(
        title = paste("ANOVA Plot:", out_variable),
        x = factor_grouping,
        y = out_variable
      )
    p
  }, error = function(e) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste("Error creating plot:", e$message), size = 6) +
      getStatfarmerTheme() +
      theme_void()
  })
}


