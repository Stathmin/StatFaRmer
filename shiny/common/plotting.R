# Plotting utilities

library(ggplot2)
library(rlang)
library(cowplot)

if (!exists(".sf_hortcontrol_trait_label", mode = "function")) {
  .hc <- file.path("shiny", "common", "hortcontrol_labels.R")
  if (!file.exists(.hc) && exists("here", mode = "function")) {
    .hc <- here::here("shiny", "common", "hortcontrol_labels.R")
  }
  if (file.exists(.hc)) source(.hc, local = FALSE)
}

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

#' Genotype count above which per-cultivar legends / violins become unreadable.
.sf_plot_high_card <- 24L

#' Keep top/bottom ends (+ checks) so ~325-genotype forests stay readable.
.sf_keep_rank_ends <- function(df, rank_col = "rank", level_col = "level",
                               keep_ends = 25L, checks = character(0)) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 1L) return(df)
  keep_ends <- max(5L, as.integer(keep_ends))
  if (nrow(df) <= 2L * keep_ends + length(checks)) return(df)
  if (!rank_col %in% names(df)) {
    df[[rank_col]] <- seq_len(nrow(df))
  }
  top <- df[df[[rank_col]] <= keep_ends, , drop = FALSE]
  max_r <- max(df[[rank_col]], na.rm = TRUE)
  bot <- df[df[[rank_col]] > (max_r - keep_ends), , drop = FALSE]
  chk <- character(0)
  if (length(checks) && level_col %in% names(df)) {
    chk <- df[as.character(df[[level_col]]) %in% as.character(checks), , drop = FALSE]
  }
  out <- rbind(top, bot, chk)
  out[!duplicated(as.character(out[[level_col]])), , drop = FALSE]
}

#' Forest / ranked-dot plot for genotype BLUPs (± SE).
#' @param genotype_table data.frame with level, blup, se, rank
#' @param checks optional check genotype names to highlight
#' @param title plot title
#' @param keep_ends if not NULL and n > 2*keep_ends, show only top/bottom ends (+ checks)
#' @return ggplot
createBlupRankPlot <- function(genotype_table, checks = NULL, title = "Genotype BLUP ranks",
                               keep_ends = NULL) {
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
  checks <- as.character(checks)
  checks <- checks[!is.na(checks) & nzchar(checks)]
  n_full <- nrow(df)
  if (!is.null(keep_ends)) {
    df <- .sf_keep_rank_ends(df, keep_ends = keep_ends, checks = checks)
  }
  df$level <- factor(df$level, levels = rev(df$level[order(df$rank, df$level)]))
  df$is_check <- df$level %in% checks
  df$lo <- df$blup - df$se
  df$hi <- df$blup + df$se
  y_cex <- if (nrow(df) > 40L) 6.5 else if (nrow(df) > 20L) 8 else 10
  sub <- NULL
  if (length(checks)) sub <- paste("Checks:", paste(checks, collapse = ", "))
  if (!is.null(keep_ends) && n_full > nrow(df)) {
    note <- sprintf("Showing top/bottom %d of %d genotypes (+ checks)", as.integer(keep_ends), n_full)
    sub <- if (is.null(sub)) note else paste(sub, "|", note)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = blup, y = level, color = is_check)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi),
      orientation = "y", width = 0.25, na.rm = TRUE, alpha = 0.7
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#2c3e50", `TRUE` = "#e67e22"),
      labels = c(`FALSE` = "genotype", `TRUE` = "check"),
      name = NULL
    ) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = sub,
      x = "BLUP",
      y = NULL
    ) +
    ggplot2::theme(
      legend.position = if (length(checks)) "bottom" else "none",
      axis.text.y = ggplot2::element_text(size = y_cex),
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
  attr(p, "sf_n_rows_plotted") <- nrow(df)
  attr(p, "sf_n_rows_full") <- n_full
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
#' @param max_show cap rows for ~325-G screens (NULL = all)
#' @param prefer_sig keep all p<0.05 first, then largest |estimate|
#' @return ggplot
createBlupDunnettPlot <- function(dunnett_df, title = "Dunnett vs check (FDR)",
                                  max_show = NULL, prefer_sig = TRUE) {
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
  n_full <- nrow(df)
  if (!is.null(max_show) && n_full > as.integer(max_show)) {
    max_show <- as.integer(max_show)
    if (isTRUE(prefer_sig)) {
      sig <- df[isTRUE(df$.sig) | (!is.na(df$.sig) & df$.sig), , drop = FALSE]
      rest <- df[is.na(df$.sig) | !df$.sig, , drop = FALSE]
      rest <- rest[order(-abs(rest$.estimate)), , drop = FALSE]
      need <- max(0L, max_show - nrow(sig))
      df <- rbind(sig, rest[seq_len(min(need, nrow(rest))), , drop = FALSE])
      if (nrow(df) > max_show) {
        df <- df[order(-abs(df$.estimate)), , drop = FALSE]
        df <- df[seq_len(max_show), , drop = FALSE]
      }
    } else {
      df <- df[order(-abs(df$.estimate)), , drop = FALSE]
      df <- df[seq_len(max_show), , drop = FALSE]
    }
  }
  df$.label <- factor(df$.label, levels = rev(df$.label[order(df$.estimate)]))
  df$lo <- df$.estimate - df$.se
  df$hi <- df$.estimate + df$.se
  ref <- attr(dunnett_df, "dunnett_ref")
  y_cex <- if (nrow(df) > 40L) 6.5 else if (nrow(df) > 20L) 8 else 10
  sub <- if (!is.null(ref) && nzchar(as.character(ref))) {
    paste("Reference check:", ref)
  } else {
    NULL
  }
  if (n_full > nrow(df)) {
    note <- sprintf("Showing %d of %d contrasts (sig + largest |Δ|)", nrow(df), n_full)
    sub <- if (is.null(sub)) note else paste(sub, "|", note)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .estimate, y = .label, color = .sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lo, xmax = hi),
      orientation = "y", width = 0.25, na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "#7f8c8d", `TRUE` = "#c0392b"),
      labels = c(`FALSE` = "n.s.", `TRUE` = "p < 0.05 (FDR)"),
      name = NULL
    ) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = sub,
      x = "Contrast estimate",
      y = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = y_cex),
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
}

#' PCA scree (proportion of variance) — ggplot, shared theme.
#' @param variance_df data.frame with component, prop_var
#' @return ggplot
createDrScreePlot <- function(variance_df, title = "Scree (plant × traits PCA)", max_show = 15L) {
  if (is.null(variance_df) || !is.data.frame(variance_df) || nrow(variance_df) < 1L) {
    return(statfarmerEmptyPlot("No PCA variance table — Calculate PCA first."))
  }
  if (!all(c("component", "prop_var") %in% names(variance_df))) {
    return(statfarmerEmptyPlot("Variance table missing component/prop_var."))
  }
  df <- variance_df
  n_show <- min(nrow(df), max(3L, as.integer(max_show)))
  df <- df[seq_len(n_show), , drop = FALSE]
  df$pc_num <- as.integer(gsub("^PC", "", as.character(df$component)))
  if (anyNA(df$pc_num)) df$pc_num <- seq_len(nrow(df))
  df$pc_num <- factor(df$pc_num, levels = df$pc_num)
  ggplot2::ggplot(df, ggplot2::aes(x = pc_num, y = prop_var, group = 1)) +
    ggplot2::geom_col(fill = "#3498db", width = 0.75, alpha = 0.9) +
    ggplot2::geom_line(color = "#2c3e50", linewidth = 0.6) +
    ggplot2::geom_point(color = "#2c3e50", size = 2) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = if (nrow(variance_df) > n_show) {
        sprintf("Showing first %d of %d components (x = PC index)", n_show, nrow(variance_df))
      } else {
        NULL
      },
      x = "PC #",
      y = "Proportion of variance"
    ) +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))
}

#' Scale loading vectors onto a score cloud (biplot-style).
#' Arrows stay inside the score cloud (with margin) so labels are not clipped.
.sf_biplot_vectors <- function(loadings_df, x_pc, y_pc, score_df,
                               traits = c("digital_biomass_mm3", "ndvi_average"),
                               scale = 0.95) {
  if (is.null(loadings_df) || !is.data.frame(loadings_df) || !"trait" %in% names(loadings_df)) {
    return(NULL)
  }
  if (!x_pc %in% names(loadings_df) || !y_pc %in% names(loadings_df)) return(NULL)
  want <- as.character(traits)
  have <- as.character(loadings_df$trait)
  # Prefer exact names; else resolve base trait → __auc then __last
  resolved <- character(0)
  for (w in want) {
    if (w %in% have) {
      resolved <- c(resolved, w)
      next
    }
    hits <- grep(paste0("^", gsub("([\\W])", "\\\\\\1", w), "__(auc|last|peak|growth_rate|mean)$"),
                 have, value = TRUE, ignore.case = TRUE)
    if (!length(hits)) {
      hits <- grep(w, have, value = TRUE, ignore.case = TRUE)
    }
    # One arrow per base trait: prefer auc, then last
    if (length(hits)) {
      pick <- hits[grepl("__auc$", hits)][1]
      if (is.na(pick)) pick <- hits[grepl("__last$", hits)][1]
      if (is.na(pick)) pick <- hits[[1L]]
      resolved <- c(resolved, pick)
    }
  }
  resolved <- unique(resolved[!is.na(resolved)])
  if (length(resolved) < 1L) return(NULL)
  ld <- loadings_df[match(resolved, have), , drop = FALSE]
  if (nrow(ld) < 1L) return(NULL)
  sx <- as.numeric(score_df[[x_pc]])
  sy <- as.numeric(score_df[[y_pc]])
  sx <- sx[is.finite(sx)]
  sy <- sy[is.finite(sy)]
  if (length(sx) < 2L || length(sy) < 2L) return(NULL)
  # Size arrows to the dense cloud (95% |score|), not outlier tails that
  # otherwise shrink arrows into a speck — or clip them on the short side.
  sc <- max(0.35, min(1.0, as.numeric(scale)))
  qx <- as.numeric(stats::quantile(abs(sx), probs = 0.95, names = FALSE, na.rm = TRUE))
  qy <- as.numeric(stats::quantile(abs(sy), probs = 0.95, names = FALSE, na.rm = TRUE))
  reach_x <- qx * sc
  reach_y <- qy * sc
  if (!is.finite(reach_x) || !is.finite(reach_y) || reach_x < 1e-9 || reach_y < 1e-9) {
    return(NULL)
  }
  vx <- as.numeric(ld[[x_pc]])
  vy <- as.numeric(ld[[y_pc]])
  need <- max(c(abs(vx) / reach_x, abs(vy) / reach_y), na.rm = TRUE)
  if (!is.finite(need) || need < 1e-9) return(NULL)
  mult <- 1 / need
  data.frame(
    trait = as.character(ld$trait),
    label = vapply(as.character(ld$trait), function(nm) {
      if (exists(".sf_hortcontrol_trait_label", mode = "function")) {
        .sf_hortcontrol_trait_label(nm)
      } else {
        nm
      }
    }, character(1)),
    x0 = 0, y0 = 0,
    x1 = vx * mult,
    y1 = vy * mult,
    stringsAsFactors = FALSE
  )
}

#' Plant PC score scatter (PC1 vs PC2 by default), colored by cultivar/treatment.
#' @param plant_df plant table with PC* columns
#' @param color_by aesthetic factor column if present
#' @param x_pc,y_pc PC column names
#' @param facet_fml optional facet formula (uses applyStatfarmerFacet)
#' @param loadings optional loadings for biplot arrows
#' @param vector_traits trait names to draw as arrows (default biomass + NDVI)
#' @param subtitle optional subtitle override / append
#' @return ggplot
createDrScoresPlot <- function(plant_df,
                               color_by = "cultivar",
                               x_pc = "PC1",
                               y_pc = "PC2",
                               facet_fml = NULL,
                               facet_ncol = 2L,
                               title = "Plant PC scores",
                               loadings = NULL,
                               vector_traits = c("digital_biomass_mm3", "ndvi_average"),
                               subtitle = NULL) {
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
  nlev_color <- if (!is.null(color_by) && color_by %in% names(df)) {
    length(unique(stats::na.omit(df[[color_by]])))
  } else {
    0L
  }
  # ~325 cultivars: a color legend is unreadable — grey cloud instead
  use_color <- nlev_color >= 1L && nlev_color <= .sf_plot_high_card
  high_card_note <- if (nlev_color > .sf_plot_high_card) {
    sprintf("%d %ss (no color legend)", nlev_color, color_by)
  } else {
    NULL
  }
  vec <- if (!identical(x_pc, y_pc)) {
    .sf_biplot_vectors(loadings, x_pc, y_pc, df, traits = vector_traits)
  } else {
    NULL
  }
  vec_note <- if (!is.null(vec) && nrow(vec) > 0L) {
    paste0("Arrows: ", paste(unique(vec$label), collapse = " + "))
  } else {
    NULL
  }
  # Keep top subtitle short; park meta (arrows / high-card) in caption so it fits
  sub_final <- if (!is.null(subtitle) && nzchar(as.character(subtitle)[[1L]])) {
    as.character(subtitle)[[1L]]
  } else {
    NULL
  }
  cap_bits <- c(vec_note, high_card_note)
  cap_bits <- cap_bits[!vapply(cap_bits, function(z) is.null(z) || !nzchar(as.character(z)), logical(1))]
  cap_final <- if (length(cap_bits)) paste(cap_bits, collapse = " · ") else NULL

  if (identical(x_pc, y_pc)) {
    # Single PC: strip / jitter vs index
    df$.idx <- seq_len(nrow(df))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .idx, y = !!rlang::sym(x_pc)))
    if (use_color) {
      p <- p + ggplot2::aes(color = !!rlang::sym(color_by))
    }
    if (use_color) {
      p <- p + ggplot2::geom_point(size = 2, alpha = 0.85)
    } else {
      p <- p + ggplot2::geom_point(size = 2, alpha = 0.45, color = "#2c3e50")
    }
    p <- p +
      getStatfarmerTheme() +
      ggplot2::labs(title = title, subtitle = sub_final, caption = cap_final,
                    x = "Plant (complete cases)", y = x_pc)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(x_pc), y = !!rlang::sym(y_pc)))
    if (use_color) {
      p <- p + ggplot2::aes(color = !!rlang::sym(color_by))
    }
    p <- p +
      ggplot2::geom_hline(yintercept = 0, color = "grey80") +
      ggplot2::geom_vline(xintercept = 0, color = "grey80")
    if (use_color) {
      p <- p + ggplot2::geom_point(size = 2.0, alpha = 0.85)
    } else {
      p <- p + ggplot2::geom_point(size = 2.0, alpha = 0.4, color = "#2c3e50")
    }
    if (!is.null(vec) && nrow(vec) > 0L) {
      sx <- df[[x_pc]]
      sy <- df[[y_pc]]
      # Full score range so every point sits inside the axes (no floaters outside panel)
      xr <- range(sx, finite = TRUE)
      yr <- range(sy, finite = TRUE)
      vec$lx <- vec$x1 * 0.62
      vec$ly <- vec$y1 * 0.62
      pad_x <- diff(xr) * 0.06
      pad_y <- diff(yr) * 0.06
      xlim <- range(c(xr, 0, vec$x1), finite = TRUE) + c(-pad_x, pad_x)
      ylim <- range(c(yr, 0, vec$y1), finite = TRUE) + c(-pad_y, pad_y)
      p <- p +
        ggplot2::geom_segment(
          data = vec,
          ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
          inherit.aes = FALSE,
          arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed"),
          color = "#c0392b", linewidth = 1.05, lineend = "round"
        ) +
        ggplot2::geom_label(
          data = vec,
          ggplot2::aes(x = lx, y = ly, label = label),
          inherit.aes = FALSE,
          color = "#c0392b", fill = "white", linewidth = 0.2,
          fontface = "bold", size = 3.3, label.padding = ggplot2::unit(0.12, "lines")
        ) +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "on", expand = FALSE)
    }
    p <- p +
      getStatfarmerTheme() +
      ggplot2::labs(title = title, subtitle = sub_final, caption = cap_final,
                    x = x_pc, y = y_pc) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(16, 16, 20, 16),
        plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5,
                                           margin = ggplot2::margin(t = 2, b = 4)),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5,
                                              margin = ggplot2::margin(b = 8)),
        plot.caption = ggplot2::element_text(size = 8.5, hjust = 0.5, color = "grey30",
                                             margin = ggplot2::margin(t = 10, b = 4))
      )
  }
  if (use_color) {
    nlev <- nlev_color
    if (nlev <= 20L && requireNamespace("viridis", quietly = TRUE)) {
      p <- p + viridis::scale_color_viridis(discrete = TRUE, end = 0.9)
    }
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))
    ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9)
      )
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.null(facet_fml)) {
    p <- applyStatfarmerFacet(p, facet_fml, ncol = facet_ncol, data = df)
  }
  p
}

#' Trait loadings as long-format bars (one panel per PC) via shared theme + facet_wrap.
#' @param loadings_df wide data.frame(trait, PC1, PC2, …)
#' @param max_pcs max components to show
#' @param top_n keep top |loading| traits per PC (avoids unreadable overlap)
#' @return ggplot
createDrLoadingsPlot <- function(loadings_df, max_pcs = 3L, top_n = 12L,
                                 title = "Trait loadings") {
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
  top_n <- max(5L, as.integer(top_n))

  long <- do.call(rbind, lapply(pc_cols, function(pc) {
    data.frame(
      trait = as.character(loadings_df$trait),
      component = pc,
      loading = as.numeric(loadings_df[[pc]]),
      stringsAsFactors = FALSE
    )
  }))
  # Per-PC top |loading| only — full 40-trait panels crush y-labels
  long <- do.call(rbind, lapply(split(long, long$component), function(chunk) {
    chunk <- chunk[order(-abs(chunk$loading)), , drop = FALSE]
    chunk <- chunk[seq_len(min(top_n, nrow(chunk))), , drop = FALSE]
    chunk
  }))
  long$component <- factor(long$component, levels = pc_cols)
  # HortControl-style axis labels (GLI bin [0, 0.15] · auc, Digital Biomass · last, …)
  lab_fun <- if (exists(".sf_hortcontrol_trait_label", mode = "function")) {
    .sf_hortcontrol_trait_label
  } else {
    identity
  }
  long$trait_lab <- vapply(as.character(long$trait), lab_fun, character(1))
  # Unique y keys so free_y facets don't collide / overlap labels
  long <- long[order(long$component, long$loading), , drop = FALSE]
  long$.ylab <- paste0(as.character(long$component), " · ", long$trait_lab)
  long$.ylab <- factor(long$.ylab, levels = unique(long$.ylab))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = loading, y = .ylab, fill = loading > 0)) +
    ggplot2::geom_col(width = 0.75, alpha = 0.9, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#3498db", `FALSE` = "#e74c3c")) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^PC[0-9]+ · ", "", x)) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Top %d |loadings| per PC (HortControl names)", top_n),
      x = "Loading",
      y = NULL
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8.5),
      strip.text = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    ) +
    ggplot2::facet_wrap(~component, ncol = 1L, scales = "free_y")
  p
}

# =============================================================================
# Trait follow-up: growth curves + violins + Tukey letters (main-plot scheme)
# =============================================================================

.sf_pick_time_col <- function(data) {
  if ("timestamp_group" %in% names(data)) return("timestamp_group")
  if ("time_numeric" %in% names(data)) return("time_numeric")
  if ("dbscan_cluster" %in% names(data)) return("dbscan_cluster")
  NULL
}

#' Growth curve: raw points (faint) + median trajectory by group.
#' Same theme / viridis language as main distPlot timeseries.
#' High genotype cardinality: grey cloud + overall median (no 325-color legend).
createTraitGrowthCurvePlot <- function(data,
                                       trait,
                                       group = "cultivar",
                                       time_col = NULL,
                                       title = NULL) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
    return(statfarmerEmptyPlot("No data for growth curve."))
  }
  if (!trait %in% names(data)) {
    return(statfarmerEmptyPlot(paste("Trait missing:", trait)))
  }
  if (is.null(time_col)) time_col <- .sf_pick_time_col(data)
  if (is.null(time_col) || !time_col %in% names(data)) {
    return(statfarmerEmptyPlot("No time column (dbscan_cluster / timestamp_group)."))
  }
  use_group <- !is.null(group) && group %in% names(data)
  df <- data
  df[[trait]] <- suppressWarnings(as.numeric(df[[trait]]))
  df <- df[is.finite(df[[trait]]), , drop = FALSE]
  if (nrow(df) < 1L) return(statfarmerEmptyPlot("No finite trait values."))
  nlev <- if (use_group) length(unique(stats::na.omit(df[[group]]))) else 0L
  high_card <- use_group && nlev > .sf_plot_high_card

  if (high_card) {
    # Quantile ribbon + overall median (325 spaghetti lines are unreadable)
    band <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(time_col))) |>
      dplyr::summarise(
        .y = stats::median(.data[[trait]], na.rm = TRUE),
        .lo = stats::quantile(.data[[trait]], 0.1, na.rm = TRUE),
        .hi = stats::quantile(.data[[trait]], 0.9, na.rm = TRUE),
        .groups = "drop"
      )
    # Sample a few genotype median trajectories for texture (not a legend)
    g_ids <- unique(as.character(df[[group]]))
    set.seed(1L)
    g_sample <- if (length(g_ids) > 12L) sample(g_ids, 12L) else g_ids
    med_g <- df[as.character(df[[group]]) %in% g_sample, , drop = FALSE] |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(time_col, group)))) |>
      dplyr::summarise(.y = stats::median(.data[[trait]], na.rm = TRUE), .groups = "drop")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]], y = .data[[trait]])) +
      ggplot2::geom_ribbon(
        data = band,
        ggplot2::aes(x = .data[[time_col]], ymin = .lo, ymax = .hi),
        inherit.aes = FALSE, fill = "#aed6f1", alpha = 0.45
      ) +
      ggplot2::geom_point(
        alpha = 0.04, size = 0.7, color = "#7f8c8d",
        position = ggplot2::position_jitter(width = 0.1, height = 0)
      ) +
      ggplot2::geom_line(
        data = med_g,
        ggplot2::aes(y = .y, group = .data[[group]]),
        color = "#7f8c8d", linewidth = 0.35, alpha = 0.45
      ) +
      ggplot2::geom_line(
        data = band, ggplot2::aes(y = .y, group = 1),
        color = "#c0392b", linewidth = 1.2
      ) +
      ggplot2::geom_point(
        data = band, ggplot2::aes(y = .y),
        color = "#c0392b", size = 2.5
      ) +
      getStatfarmerTheme() +
      ggplot2::labs(
        title = title %||% paste("Growth curve:", trait),
        subtitle = sprintf(
          "%d genotypes: 10–90%% ribbon + red median (+ %d sample trajectories)",
          nlev, length(g_sample)
        ),
        x = if (identical(time_col, "timestamp_group")) "Date/Time" else if (identical(time_col, "time_numeric")) "Time (days)" else "Time cluster",
        y = trait
      ) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = if (identical(time_col, "timestamp_group")) {
          ggplot2::element_text(angle = 45, hjust = 1, size = 8)
        } else {
          ggplot2::element_text(size = 9)
        }
      )
    return(p)
  }

  if (use_group) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[time_col]], y = .data[[trait]],
        color = .data[[group]], group = .data[[group]]
      )
    )
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]], y = .data[[trait]], group = 1))
  }
  p <- p +
    ggplot2::geom_point(alpha = 0.12, size = 1.1, position = ggplot2::position_jitter(width = 0.15, height = 0)) +
    ggplot2::stat_summary(fun = stats::median, geom = "line", linewidth = 0.95) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", size = 2.2) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title %||% paste("Growth curve:", trait),
      subtitle = "Raw points + median trajectory (StatFaRmer theme)",
      x = if (identical(time_col, "timestamp_group")) "Date/Time" else if (identical(time_col, "time_numeric")) "Time (days)" else "Time cluster",
      y = trait,
      color = group
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      axis.text.x = if (identical(time_col, "timestamp_group")) {
        ggplot2::element_text(angle = 45, hjust = 1, size = 8)
      } else {
        ggplot2::element_text(size = 9)
      }
    )
  if (use_group && requireNamespace("viridis", quietly = TRUE)) {
    if (nlev <= 20L) {
      p <- p + viridis::scale_color_viridis(discrete = TRUE, end = 0.9)
    }
  }
  p + ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 3)))
}

#' High-cardinality single-time distribution: top/bottom genotype ends (no 325 violins).
#' @param ends how many highest + lowest genotype means to show
createTraitRankStripPlot <- function(data,
                                     trait,
                                     group = "cultivar",
                                     title = NULL,
                                     subtitle = NULL,
                                     ends = 40L) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
    return(statfarmerEmptyPlot("No data for rank-strip plot."))
  }
  if (!trait %in% names(data) || !group %in% names(data)) {
    return(statfarmerEmptyPlot("Trait or grouping column missing."))
  }
  df <- data
  df[[trait]] <- suppressWarnings(as.numeric(df[[trait]]))
  df <- df[is.finite(df[[trait]]) & !is.na(df[[group]]), , drop = FALSE]
  if (nrow(df) < 1L) return(statfarmerEmptyPlot("No finite values for rank-strip."))
  ends <- max(5L, as.integer(ends))
  gmean <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group))) |>
    dplyr::summarise(.gm = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.gm)) |>
    dplyr::mutate(.rk = dplyr::row_number())
  n_g <- nrow(gmean)
  if (n_g > 2L * ends) {
    gmean <- gmean |>
      dplyr::filter(.rk <= ends | .rk > (n_g - ends)) |>
      dplyr::mutate(
        .panel = dplyr::if_else(.rk <= ends, "Highest means", "Lowest means")
      )
    panel_note <- sprintf(
      "Top %d + bottom %d of %d genotypes by mean (use BLUP/Dunnett, not Tukey letters)",
      ends, ends, n_g
    )
  } else {
    gmean$.panel <- "All genotypes"
    panel_note <- sprintf("%d genotypes (all shown)", n_g)
  }
  keep <- as.character(gmean[[group]])
  df <- df[as.character(df[[group]]) %in% keep, , drop = FALSE]
  df$.panel <- gmean$.panel[match(as.character(df[[group]]), as.character(gmean[[group]]))]
  ord <- as.character(gmean[[group]])
  df[[group]] <- factor(df[[group]], levels = rev(ord))
  y_cex <- if (nrow(gmean) > 50L) 5.5 else if (nrow(gmean) > 30L) 7 else 8

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[trait]], y = .data[[group]])) +
    ggplot2::geom_boxplot(
      orientation = "y", width = 0.55, outlier.shape = NA,
      fill = "#d6eaf8", color = "#2c3e50", alpha = 0.85
    ) +
    ggplot2::geom_jitter(
      height = 0.12, width = 0, alpha = 0.35, size = 1.0, color = "#34495e"
    ) +
    ggplot2::facet_wrap(~.panel, scales = "free_y", ncol = 1L) +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title %||% paste("Ranked single-time:", trait),
      subtitle = subtitle %||% panel_note,
      x = trait,
      y = NULL
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = y_cex),
      strip.text = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 14, 10, 8)
    )
}

#' Violin + box + optional Tukey CLD letters on group (same letter overlay as distPlot).
#' Auto-switches to createTraitRankStripPlot when genotype count is high.
createTraitViolinPlot <- function(data,
                                  trait,
                                  group = "cultivar",
                                  letters_df = NULL,
                                  title = NULL,
                                  subtitle = NULL,
                                  raw_jitter = TRUE,
                                  ends = 40L) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) < 1L) {
    return(statfarmerEmptyPlot("No data for violin plot."))
  }
  if (!trait %in% names(data) || !group %in% names(data)) {
    return(statfarmerEmptyPlot("Trait or grouping column missing."))
  }
  df <- data
  df[[trait]] <- suppressWarnings(as.numeric(df[[trait]]))
  df <- df[is.finite(df[[trait]]) & !is.na(df[[group]]), , drop = FALSE]
  if (nrow(df) < 1L) return(statfarmerEmptyPlot("No finite values for violin."))
  nlev <- length(unique(stats::na.omit(df[[group]])))
  if (nlev > .sf_plot_high_card) {
    return(createTraitRankStripPlot(
      df, trait, group = group, title = title, subtitle = subtitle, ends = ends
    ))
  }
  df[[group]] <- factor(df[[group]])

  # Join Tukey letters if provided (letters_df: group, letter)
  if (!is.null(letters_df) && is.data.frame(letters_df) && nrow(letters_df) > 0L &&
      all(c("group", "letter") %in% names(letters_df))) {
    if (exists("joinLettersToData", mode = "function")) {
      jl <- joinLettersToData(df, letters_df, group)
      df <- jl$df
      let_col <- paste0("_", group, "_letter")
      if (let_col %in% names(df)) {
        df$group_letter <- as.character(df[[let_col]])
      }
    } else {
      df$group_letter <- letters_df$letter[match(as.character(df[[group]]), as.character(letters_df$group))]
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group]], y = .data[[trait]])) +
    ggplot2::geom_violin(fill = "#aed6f1", color = "#2c3e50", alpha = 0.55, scale = "width") +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, fill = "white") +
    getStatfarmerTheme() +
    ggplot2::labs(
      title = title %||% paste("Distribution:", trait),
      subtitle = subtitle %||% (
        if ("group_letter" %in% names(df) && any(!is.na(df$group_letter))) {
          "Single-time plants + Tukey CLD (not an across-time mean)"
        } else {
          "Single-time plants + boxplots (not an across-time mean)"
        }
      ),
      x = group,
      y = trait
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 9))

  if (isTRUE(raw_jitter)) {
    p <- p + ggplot2::geom_jitter(
      width = 0.12, height = 0, alpha = 0.25, size = 1.1, color = "#34495e"
    )
  }

  if ("group_letter" %in% names(df) && any(!is.na(df$group_letter))) {
    letter_pos <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group))) |>
      dplyr::summarise(
        .y = max(.data[[trait]], na.rm = TRUE),
        group_letter = {
          vals <- group_letter[!is.na(group_letter)]
          if (length(vals)) as.character(vals[[1]]) else NA_character_
        },
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(group_letter), is.finite(.y))
    if (nrow(letter_pos) > 0L) {
      y_top <- max(letter_pos$.y, na.rm = TRUE)
      p <- p +
        ggplot2::geom_text(
          data = letter_pos,
          ggplot2::aes(x = .data[[group]], y = .y, label = group_letter),
          vjust = -0.45, size = 4, fontface = "bold", inherit.aes = FALSE
        ) +
        ggplot2::expand_limits(y = y_top * 1.1)
    }
  }
  p
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

