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
      
      # Grid and background (explicit white — avoids black SVG / blank redraws)
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
}

#' Apply facet_grid or facet_wrap from a user formula.
#' One-sided / `. ~ x` / `x ~ .` use facet_wrap so long factor strips can wrap (#14).
#' @param p ggplot
#' @param facet_fml formula or NULL
#' @param ncol wrap columns when wrapping
#' @return ggplot
applyStatfarmerFacet <- function(p, facet_fml, ncol = 2L) {
  if (is.null(facet_fml) || inherits(facet_fml, "try-error")) return(p)
  rhs <- tryCatch(facet_fml[[length(facet_fml)]], error = function(e) NULL)
  lhs <- if (length(facet_fml) >= 3) facet_fml[[2]] else NULL
  one_sided <- length(facet_fml) == 2
  lhs_dot <- !is.null(lhs) && (identical(lhs, as.name(".")) || identical(lhs, quote(.)))
  rhs_dot <- !is.null(rhs) && (identical(rhs, as.name(".")) || identical(rhs, quote(.)))
  if (one_sided || lhs_dot || rhs_dot) {
    wrap_fml <- if (one_sided) facet_fml else if (lhs_dot) stats::as.formula(paste("~", deparse(rhs))) else stats::as.formula(paste("~", deparse(lhs)))
    p + ggplot2::facet_wrap(wrap_fml, ncol = max(1L, as.integer(ncol)))
  } else {
    p + ggplot2::facet_grid(facet_fml)
  }
}

#' Create ANOVA plot
createANOVAPlot <- function(data, formula, out_variable, factor_grouping) {
  tryCatch({
    if (is.null(factor_grouping) || !nzchar(factor_grouping) || identical(factor_grouping, "none")) {
      data$.sf_x <- "all"
      x_var <- ".sf_x"
      x_lab <- "all"
    } else {
      x_var <- factor_grouping
      x_lab <- factor_grouping
    }
    p <- data %>%
      ggplot(aes(x = !!sym(x_var), y = !!sym(out_variable))) +
      geom_boxplot() +
      getStatfarmerTheme() +
      labs(
        title = paste("ANOVA Plot:", out_variable),
        x = x_lab,
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


