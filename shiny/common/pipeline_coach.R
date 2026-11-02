# Helpers for the Clippy-style analysis coach (GS → Trait DR → BLUP on PCs).

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Explain a PC from absolute loadings (human-readable).
#' @param loadings data.frame with trait + PC columns
#' @param pc column name e.g. "PC1"
#' @param top_n how many traits to list
#' @return character scalar
.sf_explain_pc <- function(loadings, pc, top_n = 3L) {
  if (is.null(loadings) || !is.data.frame(loadings) || !"trait" %in% names(loadings)) {
    return("No loadings yet — run Trait DR (Calculate PCA).")
  }
  if (!pc %in% names(loadings)) {
    return(paste0(pc, " not in loadings."))
  }
  vals <- setNames(as.numeric(loadings[[pc]]), as.character(loadings$trait))
  ord <- order(abs(vals), decreasing = TRUE)
  n <- min(as.integer(top_n), length(ord))
  top_names <- names(vals)[ord[seq_len(n)]]
  bits <- sprintf("%s (%+.2f)", top_names, vals[top_names])

  axis_hint <- "trait mixture"
  if (all(c("mean", "auc", "peak") %in% names(vals))) {
    v <- abs(vals[c("mean", "auc", "peak")])
    if ((max(v) - min(v)) < 0.08) {
      axis_hint <- "overall curve size (mean ≈ AUC ≈ peak)"
    } else if (sign(vals[["peak"]]) != sign(vals[["auc"]]) &&
               abs(vals[["peak"]]) > 0.4 && abs(vals[["auc"]]) > 0.4) {
      axis_hint <- "curve shape: peak vs integrated AUC trade-off"
    } else {
      axis_hint <- "weighted mix of curve metrics"
    }
  }
  paste0(axis_hint, " — top: ", paste(bits, collapse = "; "))
}

#' Extract Dunnett rows with p < alpha (FDR column preferred).
.sf_dunnett_significant <- function(dunnett_df, alpha = 0.05) {
  if (is.null(dunnett_df) || !is.data.frame(dunnett_df) || nrow(dunnett_df) < 1L) {
    return(if (is.data.frame(dunnett_df)) dunnett_df[0, , drop = FALSE] else data.frame())
  }
  pcol <- intersect(c("p.value", "p.adj", "adj.p.value"), names(dunnett_df))[1]
  if (is.na(pcol)) {
    return(dunnett_df[0, , drop = FALSE])
  }
  keep <- is.finite(dunnett_df[[pcol]]) & dunnett_df[[pcol]] < alpha
  dunnett_df[keep, , drop = FALSE]
}

#' Coach pipeline status from module outputs (NULL-safe).
.sf_coach_status <- function(plant_summaries = NULL,
                             dr_result = NULL,
                             blup_result = NULL,
                             alpha = 0.05) {
  steps <- c(
    gs = !is.null(plant_summaries) && is.data.frame(plant_summaries) && nrow(plant_summaries) > 0L,
    dr = !is.null(dr_result) && isTRUE(dr_result$meta$success),
    blup = !is.null(blup_result) && isTRUE(blup_result$success)
  )

  pcs <- list()
  if (isTRUE(steps[["dr"]])) {
    pcs$PC1 <- .sf_explain_pc(dr_result$loadings, "PC1")
    if ("PC2" %in% names(dr_result$loadings)) {
      pcs$PC2 <- .sf_explain_pc(dr_result$loadings, "PC2")
    }
    if (is.data.frame(dr_result$variance) && nrow(dr_result$variance) > 0L) {
      pcs$var_pc1 <- dr_result$variance$prop_var[
        dr_result$variance$component == "PC1"
      ][1]
    }
    n_tr <- length(dr_result$meta$trait_cols %||% character(0))
    pcs$n_traits <- n_tr
    if (!is.null(pcs$var_pc1) && is.finite(pcs$var_pc1)) {
      pcs$var_note <- if (!is.na(n_tr) && n_tr <= 3L) {
        sprintf(
          "PC1 ≈ %.0f%% of variance among the %d selected summaries only — not the full phenotype space. Add more traits in Trait DR.",
          100 * pcs$var_pc1, n_tr
        )
      } else {
        sprintf(
          "PC1 ≈ %.0f%% of variance among %d selected plant-level traits (still not raw time×trait wide matrix unless those columns were included).",
          100 * pcs$var_pc1, n_tr
        )
      }
    }
  }

  sig <- NULL
  ranks_head <- NULL
  trait_used <- NULL
  check_used <- NULL
  biology <- NULL
  if (isTRUE(steps[["blup"]])) {
    trait_used <- blup_result$meta$trait
    check_used <- blup_result$meta$checks
    ranks_head <- utils::head(blup_result$genotype_table, 8L)
    sig <- .sf_dunnett_significant(blup_result$meta$dunnett, alpha = alpha)

    # Prefer multi-PC biology when caller attached blup_by_pc on DR or BLUP meta
    blup_by_pc <- blup_result$meta$blup_by_pc
    if (is.null(blup_by_pc) && isTRUE(steps[["dr"]])) {
      blup_by_pc <- dr_result$meta$blup_by_pc
    }
    if (is.null(blup_by_pc) && isTRUE(steps[["blup"]])) {
      # Single-PC fallback: wrap current BLUP
      pc_name <- trait_used %||% "PC1"
      blup_by_pc <- setNames(list(blup_result), pc_name)
    }
    if (exists("buildBiologyReadout", mode = "function") && isTRUE(steps[["dr"]])) {
      biology <- tryCatch(
        buildBiologyReadout(
          loadings = dr_result$loadings,
          variance = dr_result$variance,
          blup_by_pc = blup_by_pc,
          check = check_used,
          top_k = 8L,
          overlap_k = 25L
        ),
        error = function(e) NULL
      )
    }
  }

  if (!steps[["gs"]]) {
    step <- 1L
    message <- "Hi! I'm FarmaClip — let's find which genotypes stand out on growth-curve dimensions."
    next_action <- "Open Growth Summaries → enable mean/AUC/peak → Calculate Summaries (uses left-bar times/cultivars)."
  } else if (!steps[["dr"]]) {
    step <- 2L
    message <- sprintf(
      "Growth Summaries ready (%d plants). Next we compress curve metrics into PCs.",
      nrow(plant_summaries)
    )
    next_action <- "Open Trait DR → select mean/AUC/peak → Calculate PCA. I'll explain what PC1/PC2 mean."
  } else if (!steps[["blup"]]) {
    step <- 3L
    message <- paste0(
      "PCs ready. ",
      if (!is.null(pcs$PC1)) paste0("PC1 ≈ ", pcs$PC1, ". ") else "",
      "Now rank genotypes with BLUP on a PC."
    )
    next_action <- "Open BLUP → set trait to PC1 (or PC2) → set a check cultivar → Calculate BLUP. Plots + Dunnett show who differs from the check."
  } else {
    step <- 4L
    n_sig <- if (is.data.frame(sig)) nrow(sig) else 0L
    if (n_sig > 0L) {
      message <- sprintf(
        "Done: %d genotype(s) differ from the check on %s (Dunnett FDR p < %.2f).",
        n_sig, trait_used %||% "trait", alpha
      )
    } else {
      message <- sprintf(
        "BLUP on %s finished. No Dunnett FDR hits vs check — use ranks (top BLUPs) as the ranking answer; try another check/PC or filters.",
        trait_used %||% "trait"
      )
    }
    next_action <- "Open BLUP → Curves & Tukey: pick trait + ONE time slice → Plot curves + Tukey (never across-time raw means). Also Dunnett / Answer tabs."
    if (!is.null(biology) && length(biology$bullets)) {
      message <- paste0(
        message, " Biological takeaway: ",
        paste(utils::head(biology$bullets, 4L), collapse = " · ")
      )
    }
  }

  list(
    step = step,
    steps_done = steps,
    message = message,
    next_action = next_action,
    pcs = pcs,
    sig = sig,
    ranks_head = ranks_head,
    trait_used = trait_used,
    check_used = check_used,
    biology = biology
  )
}
