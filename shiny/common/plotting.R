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
#' Missing facet vars (e.g. `treatment` on single-treatment soy) are dropped so
#' ggplot does not abort with combine_vars.
#' @param p ggplot
#' @param facet_fml formula or NULL
#' @param ncol wrap columns when wrapping
#' @param data optional data.frame used to validate facet vars
#' @return ggplot
applyStatfarmerFacet <- function(p, facet_fml, ncol = 2L, data = NULL) {
  if (is.null(facet_fml) || inherits(facet_fml, "try-error")) return(p)

  if (!is.null(data) && is.data.frame(data)) {
    vars <- setdiff(all.vars(facet_fml), ".")
    missing <- setdiff(vars, names(data))
    if (length(missing) > 0L) {
      present <- intersect(vars, names(data))
      if (length(present) == 0L) return(p)
      facet_fml <- stats::as.formula(paste("~", paste(present, collapse = " + ")))
    }
  }

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
    statfarmerEmptyPlot(paste("Error creating plot:", e$message))
  })
}

#' Empty / status ggplot using the shared theme (sidecars + main plot).
#' @param message Character label
#' @return ggplot
statfarmerEmptyPlot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5,
      label = paste(strwrap(as.character(message), width = 60), collapse = "\n"),
      size = 4.5
    ) +
    getStatfarmerTheme() +
    ggplot2::theme_void()
}

# =============================================================================
# BLUP / Trait-DR sidecar plots (same theme + viridis as main distPlot)
# =============================================================================

#' Forest / ranked-dot plot for genotype BLUPs (± SE).
#' @param genotype_table data.frame with level, blup, se, rank
#' @param checks optional check genotype names to highlight
#' @param title plot title
#' @return ggplot
createBlupRankPlot <- function(genotype_table, checks = NULL, title = "Genotype BLUP ranks") {
  if (is.null(genotype_table) || !is.data.frame(genotype_table) || nrow(genotype_table) < 1L) {
    return(statfarmerEmptyPlot("No genotype BLUP ranks to plot."))
  }
  need <- c("level", "blup")
  if (!all(need %in% names(genotype_table))) {
    return(statfarmerEmptyPlot("Genotype table missing level/blup columns."))
  }
  df <- genotype_table
  if (!"se" %in% names(df)) df$se <- NA_real_
  if (!"rank" %in% names(df)) {
    df$rank <- as.integer(rank(-df$blup, ties.method = "min"))
  }
  df$level <- factor(df$level, levels = rev(df$level[order(df$rank, df$level)]))
  checks <- as.character(checks)
  checks <- checks[!is.na(checks) & nzchar(checks)]
  df$is_check <- df$level %in% checks
  df$lo <- df$blup - df$se
  df$hi <- df$blup + df$se

  p <- ggplot2::ggplot(df, ggplot2::aes(x = blup, y = level, color = is_check)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi),
      orientation = "y", width = 0.25, na.rm = TRUE, alpha = 0.7
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#2c3e50", `TRUE` = "#e67e22"),
      labels = c(`FALSE` = "genotype", `TRUE` = "check"),
      name = NULL
    ) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = if (length(checks)) paste("Checks:", paste(checks, collapse = ", ")) else NULL,
      x = "BLUP",
      y = NULL
    ) +
    ggplot2::theme(legend.position = if (length(checks)) "bottom" else "none")
  p
}

#' Treatment fixed-effect EMMeans (± SE) using shared theme.
#' @param treatment_table data.frame with level, estimate, se
#' @param title plot title
#' @return ggplot
createBlupTreatmentPlot <- function(treatment_table, title = "Treatment effects") {
  if (is.null(treatment_table) || !is.data.frame(treatment_table) || nrow(treatment_table) < 1L) {
    return(statfarmerEmptyPlot("No treatment estimates (need ≥2 treatment levels)."))
  }
  if (!all(c("level", "estimate") %in% names(treatment_table))) {
    return(statfarmerEmptyPlot("Treatment table missing level/estimate."))
  }
  df <- treatment_table
  if (!"se" %in% names(df)) df$se <- NA_real_
  if (!"rank" %in% names(df)) {
    df$rank <- as.integer(rank(-df$estimate, ties.method = "min"))
  }
  df$level <- factor(df$level, levels = df$level[order(df$rank, df$level)])
  df$lo <- df$estimate - df$se
  df$hi <- df$estimate + df$se

  ggplot2::ggplot(df, ggplot2::aes(x = level, y = estimate)) +
    ggplot2::geom_col(fill = "#3498db", width = 0.7, alpha = 0.85) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2, na.rm = TRUE) +
    getStatfarmerTheme() +
    ggplot2::labs(title = title, x = "Treatment", y = "EMMean") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Dunnett / FDR posthoc contrasts (estimate ± SE; color by significance).
#' @param dunnett_df data.frame from emmeans contrast (or NULL / message list)
#' @param title plot title
#' @return ggplot
createBlupDunnettPlot <- function(dunnett_df, title = "Dunnett vs check (FDR)") {
  if (is.null(dunnett_df)) {
    return(statfarmerEmptyPlot("No Dunnett results — set check genotypes and re-run BLUP."))
  }
  if (is.list(dunnett_df) && !is.data.frame(dunnett_df)) {
    msg <- dunnett_df$message
    if (is.null(msg) || !nzchar(as.character(msg))) msg <- "Dunnett unavailable."
    return(statfarmerEmptyPlot(msg))
  }
  if (!is.data.frame(dunnett_df) || nrow(dunnett_df) < 1L) {
    return(statfarmerEmptyPlot("Dunnett table empty."))
  }
  df <- dunnett_df
  # emmeans contrast columns vary slightly by version
  est_col <- intersect(c("estimate", "Estimate"), names(df))[1]
  se_col <- intersect(c("SE", "se", "std.error"), names(df))[1]
  p_col <- intersect(c("p.value", "p.adj", "adj.p.value"), names(df))[1]
  lab_col <- intersect(c("contrast", "Contrasts", "level"), names(df))[1]
  if (is.na(est_col) || is.na(lab_col)) {
    return(statfarmerEmptyPlot("Dunnett table missing contrast/estimate columns."))
  }
  df$.label <- as.character(df[[lab_col]])
  df$.estimate <- as.numeric(df[[est_col]])
  df$.se <- if (!is.na(se_col)) as.numeric(df[[se_col]]) else NA_real_
  df$.p <- if (!is.na(p_col)) as.numeric(df[[p_col]]) else NA_real_
  df$.sig <- !is.na(df$.p) & df$.p < 0.05
  df$.label <- factor(df$.label, levels = rev(df$.label[order(df$.estimate)]))
  df$lo <- df$.estimate - df$.se
  df$hi <- df$.estimate + df$.se
  ref <- attr(dunnett_df, "dunnett_ref")

  ggplot2::ggplot(df, ggplot2::aes(x = .estimate, y = .label, color = .sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi),
      orientation = "y", width = 0.25, na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#7f8c8d", `TRUE` = "#c0392b"),
      labels = c(`FALSE` = "n.s.", `TRUE` = "p < 0.05 (FDR)"),
      name = NULL
    ) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = if (!is.null(ref) && nzchar(as.character(ref))) {
        paste("Reference check:", ref)
      } else {
        NULL
      },
      x = "Contrast estimate",
      y = NULL
    )
}

#' PCA scree (proportion of variance) — ggplot, shared theme.
#' @param variance_df data.frame with component, prop_var
#' @return ggplot
createDrScreePlot <- function(variance_df, title = "Scree (plant × traits PCA)") {
  if (is.null(variance_df) || !is.data.frame(variance_df) || nrow(variance_df) < 1L) {
    return(statfarmerEmptyPlot("No PCA variance table — Calculate PCA first."))
  }
  if (!all(c("component", "prop_var") %in% names(variance_df))) {
    return(statfarmerEmptyPlot("Variance table missing component/prop_var."))
  }
  df <- variance_df
  df$component <- factor(df$component, levels = df$component)
  ggplot2::ggplot(df, ggplot2::aes(x = component, y = prop_var)) +
    ggplot2::geom_col(fill = "#3498db", width = 0.75, alpha = 0.9) +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "#2c3e50", linewidth = 0.6) +
    ggplot2::geom_point(color = "#2c3e50", size = 2) +
    getStatfarmerTheme() +
    ggplot2::labs(title = title, x = "Component", y = "Proportion of variance") +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x))
}

#' Plant PC score scatter (PC1 vs PC2 by default), colored by cultivar/treatment.
#' @param plant_df plant table with PC* columns
#' @param color_by aesthetic factor column if present
#' @param x_pc,y_pc PC column names
#' @param facet_fml optional facet formula (uses applyStatfarmerFacet)
#' @return ggplot
createDrScoresPlot <- function(plant_df,
                               color_by = "cultivar",
                               x_pc = "PC1",
                               y_pc = "PC2",
                               facet_fml = NULL,
                               facet_ncol = 2L,
                               title = "Plant PC scores") {
  if (is.null(plant_df) || !is.data.frame(plant_df) || nrow(plant_df) < 1L) {
    return(statfarmerEmptyPlot("No plant scores — Calculate PCA first."))
  }
  pcs <- grep("^PC[0-9]+$", names(plant_df), value = TRUE)
  if (length(pcs) < 1L) {
    return(statfarmerEmptyPlot("No PC columns on plant table."))
  }
  if (!x_pc %in% names(plant_df)) x_pc <- pcs[[1L]]
  if (!y_pc %in% names(plant_df)) {
    y_pc <- if (length(pcs) >= 2L) pcs[[2L]] else pcs[[1L]]
  }
  df <- plant_df
  ok <- is.finite(df[[x_pc]]) & is.finite(df[[y_pc]])
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 1L) {
    return(statfarmerEmptyPlot("All PC scores are NA (incomplete trait rows)."))
  }
  use_color <- !is.null(color_by) && color_by %in% names(df) &&
    length(unique(stats::na.omit(df[[color_by]]))) >= 1L

  if (identical(x_pc, y_pc)) {
    # Single PC: strip / jitter vs index
    df$.idx <- seq_len(nrow(df))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .idx, y = !!rlang::sym(x_pc)))
    if (use_color) {
      p <- p + ggplot2::aes(color = !!rlang::sym(color_by))
    }
    p <- p +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      getStatfarmerTheme() +
      ggplot2::labs(title = title, x = "Plant (complete cases)", y = x_pc)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(x_pc), y = !!rlang::sym(y_pc)))
    if (use_color) {
      p <- p + ggplot2::aes(color = !!rlang::sym(color_by))
    }
    p <- p +
      ggplot2::geom_hline(yintercept = 0, color = "grey80") +
      ggplot2::geom_vline(xintercept = 0, color = "grey80") +
      ggplot2::geom_point(size = 2.2, alpha = 0.85) +
      getStatfarmerTheme() +
      ggplot2::labs(title = title, x = x_pc, y = y_pc)
  }
  if (use_color) {
    nlev <- length(unique(stats::na.omit(df[[color_by]])))
    if (nlev <= 12L && requireNamespace("viridis", quietly = TRUE)) {
      p <- p + viridis::scale_color_viridis(discrete = TRUE, end = 0.9)
    }
  }
  if (!is.null(facet_fml)) {
    p <- applyStatfarmerFacet(p, facet_fml, ncol = facet_ncol, data = df)
  }
  p
}

#' Trait loadings as long-format bars (one panel per PC) via shared theme + facet_wrap.
#' @param loadings_df wide data.frame(trait, PC1, PC2, …)
#' @param max_pcs max components to show
#' @return ggplot
createDrLoadingsPlot <- function(loadings_df, max_pcs = 4L, title = "Trait loadings") {
  if (is.null(loadings_df) || !is.data.frame(loadings_df) || nrow(loadings_df) < 1L) {
    return(statfarmerEmptyPlot("No loadings — Calculate PCA first."))
  }
  if (!"trait" %in% names(loadings_df)) {
    return(statfarmerEmptyPlot("Loadings table missing trait column."))
  }
  pc_cols <- grep("^PC[0-9]+$", names(loadings_df), value = TRUE)
  if (length(pc_cols) < 1L) {
    return(statfarmerEmptyPlot("Loadings table has no PC columns."))
  }
  pc_cols <- pc_cols[seq_len(min(length(pc_cols), as.integer(max_pcs)))]
  long <- do.call(rbind, lapply(pc_cols, function(pc) {
    data.frame(
      trait = as.character(loadings_df$trait),
      component = pc,
      loading = as.numeric(loadings_df[[pc]]),
      stringsAsFactors = FALSE
    )
  }))
  long$component <- factor(long$component, levels = pc_cols)
  long$trait <- factor(long$trait, levels = rev(sort(unique(as.character(long$trait)))))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = loading, y = trait, fill = loading > 0)) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#3498db", `FALSE` = "#e74c3c")) +
    getStatfarmerTheme() +
    ggplot2::labs(title = title, x = "Loading", y = NULL)
  applyStatfarmerFacet(p, ~component, ncol = min(2L, length(pc_cols)), data = long)
}

