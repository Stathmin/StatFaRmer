# Genotype BLUP ranking via lme4 (approach B stage 2).

.sf_blup_empty_geno <- function() {
  data.frame(
    level = character(0),
    blup = numeric(0),
    se = numeric(0),
    rank = integer(0),
    stringsAsFactors = FALSE
  )
}

.sf_blup_empty_trt <- function() {
  data.frame(
    level = character(0),
    estimate = numeric(0),
    se = numeric(0),
    rank = integer(0),
    stringsAsFactors = FALSE
  )
}

.sf_blup_plant_id <- function(data) {
  if ("unit" %in% names(data)) {
    return("unit")
  }
  if ("V.T.R" %in% names(data)) {
    return("V.T.R")
  }
  NULL
}

.sf_blup_time_col <- function(data) {
  if ("dbscan_cluster" %in% names(data)) {
    return("dbscan_cluster")
  }
  if ("timestamp_group" %in% names(data)) {
    return("timestamp_group")
  }
  NULL
}

#' Max distinct time values per plant (1 = already collapsed / single scan).
.sf_blup_n_time_per_plant <- function(data) {
  id_col <- .sf_blup_plant_id(data)
  time_col <- .sf_blup_time_col(data)
  if (is.null(id_col) || is.null(time_col)) {
    return(1L)
  }
  counts <- tapply(
    data[[time_col]],
    data[[id_col]],
    function(x) length(unique(x[!is.na(x)]))
  )
  if (length(counts) < 1L) {
    return(1L)
  }
  as.integer(max(counts, na.rm = TRUE))
}

.sf_blup_fail <- function(reason, n_time_per_plant = NA_integer_, meta_extra = list()) {
  meta <- c(
    list(
      success = FALSE,
      reason = reason,
      n_time_per_plant = n_time_per_plant,
      include_table_re = NA,
      table_re_note = NA_character_,
      spatial_applied = FALSE,
      treatment_n_levels = NA_integer_
    ),
    meta_extra
  )
  list(
    success = FALSE,
    reason = reason,
    genotype_table = .sf_blup_empty_geno(),
    treatment_table = .sf_blup_empty_trt(),
    varcomp = data.frame(
      grp = character(0),
      vcov = numeric(0),
      stringsAsFactors = FALSE
    ),
    formula = NA_character_,
    meta = meta,
    h2 = NA_real_
  )
}

.sf_blup_ranef_se <- function(fit, group) {
  rr <- lme4::ranef(fit, condVar = TRUE)[[group]]
  if (is.null(rr)) {
    return(rep(NA_real_, 0L))
  }
  pv <- attr(rr, "postVar")
  if (is.null(pv)) {
    return(rep(NA_real_, nrow(rr)))
  }
  # postVar is 1 x 1 x n for scalar RE
  if (length(dim(pv)) == 3L) {
    return(sqrt(as.numeric(pv[1L, 1L, ])))
  }
  rep(NA_real_, nrow(rr))
}

.sf_blup_varcomp_df <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  data.frame(
    grp = as.character(vc$grp),
    vcov = as.numeric(vc$vcov),
    stringsAsFactors = FALSE
  )
}

.sf_blup_h2 <- function(varcomp, genotype_var) {
  vg <- varcomp$vcov[varcomp$grp == genotype_var]
  ve <- varcomp$vcov[varcomp$grp == "Residual"]
  if (length(vg) != 1L || length(ve) != 1L) {
    return(NA_real_)
  }
  if (!is.finite(vg) || !is.finite(ve) || (vg + ve) <= 0) {
    return(NA_real_)
  }
  as.numeric(vg / (vg + ve))
}

#' Dunnett vs checks on a secondary fixed-genotype model (not the BLUP RE fit).
#'
#' Ranking uses `(1|G)` BLUPs; planned tests need genotype as FE so emmeans
#' Dunnett is well-defined. Fits `.y ~ .g [+ .trt] [+ (1|.rep)] [+ (1|.tab)]`.
.sf_blup_dunnett <- function(model_df,
                             checks,
                             genotype_var,
                             has_trt = FALSE,
                             has_rep = FALSE,
                             include_table_re = FALSE) {
  if (is.null(checks) || length(checks) < 1L) {
    return(NULL)
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    return(list(message = "emmeans not available; Dunnett skipped"))
  }
  checks <- as.character(checks)
  lev <- levels(model_df$.g)
  present <- intersect(checks, lev)
  if (length(present) < 1L) {
    return(list(message = "checks not found among genotype levels; Dunnett skipped"))
  }
  ref <- present[[1L]]

  rhs <- ".g"
  if (isTRUE(has_trt) && ".trt" %in% names(model_df)) {
    rhs <- paste(rhs, "+ .trt")
  }
  use_lmer <- FALSE
  if (isTRUE(has_rep) && ".rep" %in% names(model_df) &&
    nlevels(model_df$.rep) >= 2L) {
    rhs <- paste(rhs, "+ (1|.rep)")
    use_lmer <- TRUE
  }
  if (isTRUE(include_table_re) && ".tab" %in% names(model_df) &&
    nlevels(model_df$.tab) >= 2L) {
    rhs <- paste(rhs, "+ (1|.tab)")
    use_lmer <- TRUE
  }
  frm <- stats::as.formula(paste(".y ~", rhs))

  fit_fe <- tryCatch(
    {
      if (use_lmer) {
        ctrl <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
        suppressWarnings(
          lme4::lmer(frm, data = model_df, REML = TRUE, control = ctrl)
        )
      } else {
        stats::lm(frm, data = model_df)
      }
    },
    error = function(e) e
  )
  if (inherits(fit_fe, "error")) {
    return(list(
      message = paste0("Dunnett FE fit failed: ", conditionMessage(fit_fe))
    ))
  }

  emm <- tryCatch(
    {
      if (inherits(fit_fe, "lmerMod")) {
        # Prefer Satterthwaite when lmerTest is present; else asymptotic
        # (avoids tiny residual dfs that make FDR contrasts unusable).
        df_mode <- if (requireNamespace("lmerTest", quietly = TRUE)) {
          "satterthwaite"
        } else {
          "asymptotic"
        }
        emmeans::emmeans(fit_fe, specs = ".g", lmer.df = df_mode)
      } else {
        emmeans::emmeans(fit_fe, specs = ".g")
      }
    },
    error = function(e) e
  )
  if (inherits(emm, "error")) {
    return(list(message = paste0("Dunnett emmeans failed: ", conditionMessage(emm))))
  }

  # emmeans dunnett: ref can be level name or index
  ref_idx <- match(ref, as.character(emm@grid$.g))
  ctr <- tryCatch(
    emmeans::contrast(emm, method = "dunnett", ref = ref_idx, adjust = "fdr"),
    error = function(e) e
  )
  if (inherits(ctr, "error")) {
    return(list(message = paste0("Dunnett contrast failed: ", conditionMessage(ctr))))
  }
  out <- as.data.frame(ctr)
  attr(out, "dunnett_ref") <- ref
  attr(out, "dunnett_model") <- "fixed_genotype"
  attr(out, "dunnett_formula") <- paste(deparse(frm, width.cutoff = 500L), collapse = " ")
  out
}

#' Fit genotype BLUP LMM (optional spatial residualize first).
#'
#' Blocks multi-time uncollapsed plant × time panels. Optionally residualizes
#' via [residualizeSpatial], then fits `lme4::lmer` with `(1|genotype)`,
#' treatment FE when ≥2 levels, `(1|repeat)`, and `(1|table_no)` (kept after
#' spatial unless singular — documented in `meta`).
#'
#' @param data Plant-level (collapsed) data.frame.
#' @param trait Response column name.
#' @param genotype_var Genotype factor column (default `"cultivar"`).
#' @param treatment_var Treatment column (default `"treatment"`).
#' @param repeat_var Repeat / replicate column (default `"repeat"`).
#' @param table_var Table / batch column (default `"table_no"`).
#' @param apply_spatial If `TRUE`, call `residualizeSpatial` when available.
#' @param checks Optional check genotype levels for Dunnett; `NULL` → top-k only.
#' @param top_k Max ranks retained in genotype table (default `20`).
#' @return List with `genotype_table`, `treatment_table`, `varcomp`, `formula`,
#'   `meta`, `h2`, plus `success` / `reason` on failure.
fitBlupModels <- function(data,
                          trait,
                          genotype_var = "cultivar",
                          treatment_var = "treatment",
                          repeat_var = "repeat",
                          table_var = "table_no",
                          apply_spatial = TRUE,
                          checks = NULL,
                          top_k = 20) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  if (missing(trait) || length(trait) != 1L || !is.character(trait)) {
    stop("trait must be a single column name", call. = FALSE)
  }
  if (!trait %in% names(data)) {
    stop("Can't fit BLUP: missing trait column `", trait, "`", call. = FALSE)
  }
  if (!genotype_var %in% names(data)) {
    stop("Can't fit BLUP: missing genotype column `", genotype_var, "`", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(.sf_blup_fail("lme4 not available"))
  }

  n_time <- .sf_blup_n_time_per_plant(data)
  if (isTRUE(n_time > 1L)) {
    return(.sf_blup_fail(
      reason = paste0(
        "Multiple time rows per plant (max ", n_time,
        "); collapse via Growth Summaries over selected times before BLUP"
      ),
      n_time_per_plant = n_time
    ))
  }

  work <- data
  spatial_meta <- NULL
  y_col <- trait
  if (isTRUE(apply_spatial) && exists("residualizeSpatial", mode = "function")) {
    spat <- residualizeSpatial(work, y_col = trait, table_col = table_var)
    work <- spat$data
    spatial_meta <- spat$meta
    # Use residual response only when spatial stage actually applied.
    if (isTRUE(spatial_meta$spatial_applied)) {
      y_col <- spatial_meta$resid_col
    }
  }

  # Internal safe names for formula (repeat is reserved; keep originals mapped).
  model_df <- data.frame(
    .y = as.numeric(work[[y_col]]),
    .g = factor(work[[genotype_var]]),
    stringsAsFactors = FALSE
  )
  has_trt <- treatment_var %in% names(work)
  n_trt <- if (has_trt) {
    length(unique(stats::na.omit(work[[treatment_var]])))
  } else {
    0L
  }
  if (has_trt && n_trt >= 2L) {
    model_df$.trt <- factor(work[[treatment_var]])
  }
  has_rep <- repeat_var %in% names(work)
  n_rep <- if (has_rep) {
    length(unique(stats::na.omit(work[[repeat_var]])))
  } else {
    0L
  }
  if (has_rep && n_rep >= 2L) {
    model_df$.rep <- factor(work[[repeat_var]])
  }
  has_tab <- table_var %in% names(work)
  n_tab <- if (has_tab) {
    length(unique(stats::na.omit(work[[table_var]])))
  } else {
    0L
  }
  include_table_re <- has_tab && n_tab >= 2L
  table_re_note <- if (include_table_re) {
    if (isTRUE(spatial_meta$spatial_applied)) {
      paste0(
        "Keeping (1|", table_var, ") after spatial residualize ",
        "(batch RE; spatial absorbs within-table X/Y trend, not run block)."
      )
    } else {
      paste0("Including (1|", table_var, ") as run/batch random effect.")
    }
  } else {
    paste0("Omitting (1|", table_var, "): fewer than 2 table levels.")
  }
  if (include_table_re) {
    model_df$.tab <- factor(work[[table_var]])
  }

  ok <- stats::complete.cases(model_df)
  model_df <- model_df[ok, , drop = FALSE]
  if (nrow(model_df) < 3L || nlevels(model_df$.g) < 2L) {
    return(.sf_blup_fail(
      reason = "Not enough complete rows or genotype levels for BLUP",
      n_time_per_plant = n_time,
      meta_extra = list(
        spatial_applied = isTRUE(spatial_meta$spatial_applied),
        treatment_n_levels = as.integer(n_trt)
      )
    ))
  }

  rhs <- "(1|.g)"
  if (has_trt && n_trt >= 2L) {
    rhs <- paste(".trt +", rhs)
  }
  if (has_rep && n_rep >= 2L) {
    rhs <- paste(rhs, "+ (1|.rep)")
  }
  if (include_table_re) {
    rhs <- paste(rhs, "+ (1|.tab)")
  }
  frm <- stats::as.formula(paste(".y ~", rhs))
  formula_chr <- paste(deparse(frm, width.cutoff = 500L), collapse = " ")
  # Public-facing formula with original column names
  public_rhs <- paste0("(1|", genotype_var, ")")
  if (has_trt && n_trt >= 2L) {
    public_rhs <- paste(treatment_var, "+", public_rhs)
  }
  if (has_rep && n_rep >= 2L) {
    public_rhs <- paste0(public_rhs, " + (1|", repeat_var, ")")
  }
  if (include_table_re) {
    public_rhs <- paste0(public_rhs, " + (1|", table_var, ")")
  }
  public_formula <- paste(y_col, "~", public_rhs)

  ctrl <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  fit <- tryCatch(
    suppressWarnings(
      lme4::lmer(frm, data = model_df, REML = TRUE, control = ctrl)
    ),
    error = function(e) e
  )

  # If table RE is singular after spatial, drop it and refit once.
  dropped_table_for_singular <- FALSE
  if (!inherits(fit, "error") && include_table_re &&
    isTRUE(lme4::isSingular(fit, tol = 1e-4)) &&
    isTRUE(spatial_meta$spatial_applied)) {
    rhs2 <- gsub(" \\+ \\(1\\|\\.tab\\)", "", rhs, perl = TRUE)
    frm2 <- stats::as.formula(paste(".y ~", rhs2))
    fit2 <- tryCatch(
      suppressWarnings(
        lme4::lmer(frm2, data = model_df, REML = TRUE, control = ctrl)
      ),
      error = function(e) e
    )
    if (!inherits(fit2, "error")) {
      fit <- fit2
      frm <- frm2
      include_table_re <- FALSE
      dropped_table_for_singular <- TRUE
      table_re_note <- paste0(
        "Dropped (1|", table_var, ") after spatial: singular fit; ",
        "table trend already residualized."
      )
      public_formula <- sub(
        paste0(" \\+ \\(1\\|", table_var, "\\)"),
        "",
        public_formula
      )
      formula_chr <- paste(deparse(frm, width.cutoff = 500L), collapse = " ")
    }
  }

  if (inherits(fit, "error")) {
    return(.sf_blup_fail(
      reason = paste0("lmer failed: ", conditionMessage(fit)),
      n_time_per_plant = n_time,
      meta_extra = list(
        spatial_applied = isTRUE(spatial_meta$spatial_applied),
        treatment_n_levels = as.integer(n_trt),
        include_table_re = include_table_re,
        table_re_note = table_re_note,
        formula = public_formula
      )
    ))
  }

  rr <- lme4::ranef(fit, condVar = TRUE)$.g
  blups <- as.numeric(rr[, 1L])
  levels_g <- rownames(rr)
  ses <- .sf_blup_ranef_se(fit, ".g")
  # Rank 1 = highest BLUP
  rk <- rank(-blups, ties.method = "min")
  geno_tbl <- data.frame(
    level = as.character(levels_g),
    blup = blups,
    se = ses,
    rank = as.integer(rk),
    stringsAsFactors = FALSE
  )
  geno_tbl <- geno_tbl[order(geno_tbl$rank, geno_tbl$level), , drop = FALSE]
  rownames(geno_tbl) <- NULL
  top_k <- as.integer(top_k)
  if (!is.na(top_k) && top_k > 0L && nrow(geno_tbl) > top_k) {
    geno_tbl <- geno_tbl[seq_len(top_k), , drop = FALSE]
    rownames(geno_tbl) <- NULL
  }

  treatment_message <- NULL
  trt_tbl <- .sf_blup_empty_trt()
  if (!has_trt || n_trt < 2L) {
    treatment_message <- paste0(
      "Treatment panel empty: fewer than 2 levels in `", treatment_var, "`."
    )
  } else if (requireNamespace("emmeans", quietly = TRUE)) {
    df_mode_trt <- if (requireNamespace("lmerTest", quietly = TRUE)) {
      "satterthwaite"
    } else {
      "asymptotic"
    }
    emm_t <- tryCatch(
      emmeans::emmeans(fit, specs = ".trt", lmer.df = df_mode_trt),
      error = function(e) e
    )
    if (!inherits(emm_t, "error")) {
      edf <- as.data.frame(emm_t)
      trt_tbl <- data.frame(
        level = as.character(edf$.trt),
        estimate = as.numeric(edf$emmean),
        se = as.numeric(edf$SE),
        rank = as.integer(rank(-edf$emmean, ties.method = "min")),
        stringsAsFactors = FALSE
      )
      trt_tbl <- trt_tbl[order(trt_tbl$rank, trt_tbl$level), , drop = FALSE]
      rownames(trt_tbl) <- NULL
    } else {
      treatment_message <- paste0(
        "Treatment emmeans failed: ", conditionMessage(emm_t)
      )
    }
  } else {
    treatment_message <- "emmeans not available; treatment table empty."
  }

  varcomp <- .sf_blup_varcomp_df(fit)
  # Map internal grp names back for readability
  varcomp$grp[varcomp$grp == ".g"] <- genotype_var
  varcomp$grp[varcomp$grp == ".rep"] <- repeat_var
  varcomp$grp[varcomp$grp == ".tab"] <- table_var

  h2 <- .sf_blup_h2(varcomp, genotype_var)

  dunnett <- .sf_blup_dunnett(
    model_df = model_df,
    checks = checks,
    genotype_var = genotype_var,
    has_trt = has_trt && n_trt >= 2L,
    has_rep = has_rep && n_rep >= 2L,
    include_table_re = include_table_re
  )
  dunnett_message <- NULL
  if (is.list(dunnett) && !is.data.frame(dunnett) && !is.null(dunnett$message)) {
    dunnett_message <- dunnett$message
  }
  # Relabel contrast genotype column if present
  if (is.data.frame(dunnett) && ".g" %in% names(dunnett)) {
    names(dunnett)[names(dunnett) == ".g"] <- genotype_var
  }

  meta <- list(
    success = TRUE,
    n_time_per_plant = n_time,
    n_genotypes = nlevels(model_df$.g),
    n_obs = nrow(model_df),
    trait = trait,
    y_col = y_col,
    genotype_var = genotype_var,
    treatment_var = treatment_var,
    treatment_n_levels = as.integer(n_trt),
    treatment_message = treatment_message,
    include_table_re = include_table_re,
    table_re_note = table_re_note,
    dropped_table_for_singular = dropped_table_for_singular,
    spatial_applied = isTRUE(spatial_meta$spatial_applied),
    spatial_meta = spatial_meta,
    singular = isTRUE(lme4::isSingular(fit, tol = 1e-4)),
    converged = isTRUE(fit@optinfo$conv$opt == 0),
    top_k = top_k,
    checks = checks,
    dunnett = dunnett,
    dunnett_message = dunnett_message,
    message = dunnett_message
  )

  list(
    success = TRUE,
    reason = NULL,
    genotype_table = geno_tbl,
    treatment_table = trt_tbl,
    varcomp = varcomp,
    formula = public_formula,
    meta = meta,
    h2 = h2
  )
}
