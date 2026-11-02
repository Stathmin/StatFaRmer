# Plant × traits PCA (Trait DR). Rows = plants; columns = phenotypes only.
#
# Temporal HTP: phenotype columns for multi-trait DR must be *named growth
# metrics* (trait__auc / __peak / __growth_rate / __last, or single-trait GS
# columns mean/auc/peak), never silent across-time means of raw scans.
# Use collapseManyTraitsBySelectedTimes() + .sf_dr_assert_named_temporal_metrics().

.sf_dr_design_cols <- c(
  "unit", "V.T.R", "cultivar", "treatment", "repeat",
  "T:X:Y", "table_no", "pos_x", "pos_y",
  "dbscan_cluster", "timestamp_group", "n_obs", "n_times"
)

.sf_dr_pc_names <- function(names_vec) {
  grep("^PC[0-9]+$", names_vec, value = TRUE)
}

#' Eligible numeric phenotype columns (excludes design IDs and existing PCs).
#' @param plant_df Plant-level data frame
#' @return Character vector of column names
.sf_dr_eligible_traits <- function(plant_df) {
  if (is.null(plant_df) || !is.data.frame(plant_df) || ncol(plant_df) < 1L) {
    return(character(0))
  }
  num <- names(plant_df)[vapply(plant_df, is.numeric, logical(1))]
  setdiff(num, c(.sf_dr_design_cols, .sf_dr_pc_names(names(plant_df))))
}

#' Fingerprint plant table identity (ids + eligible trait names + trait values).
#' Used to detect Growth Summaries refresh so stale PC frames are not preferred.
#' @param df Plant-level data frame (with or without PC columns)
#' @return Character scalar (NA if empty)
.sf_dr_plant_fingerprint <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 1L) {
    return(NA_character_)
  }
  id_col <- if ("unit" %in% names(df)) {
    "unit"
  } else if ("V.T.R" %in% names(df)) {
    "V.T.R"
  } else {
    NULL
  }
  ids <- if (!is.null(id_col)) {
    as.character(df[[id_col]])
  } else {
    as.character(seq_len(nrow(df)))
  }
  traits <- sort(.sf_dr_eligible_traits(df))
  trait_payload <- if (length(traits) < 1L) {
    ""
  } else {
    paste(
      vapply(
        traits,
        function(nm) {
          paste(format(df[[nm]], digits = 8, scientific = TRUE), collapse = ",")
        },
        character(1)
      ),
      collapse = "|"
    )
  }
  paste(
    paste(ids, collapse = "\1"),
    paste(traits, collapse = "\1"),
    trait_payload,
    sep = "\2"
  )
}

#' TRUE when PC-augmented frame still matches current plant summaries.
#' @param plant_summaries_df Current Growth Summaries plant table
#' @param plant_with_pcs Cached PCA plant table (may include PC*)
#' @return Logical
.sf_dr_pcs_aligned <- function(plant_summaries_df, plant_with_pcs) {
  if (is.null(plant_with_pcs) || !is.data.frame(plant_with_pcs) ||
      nrow(plant_with_pcs) < 1L) {
    return(FALSE)
  }
  if (is.null(plant_summaries_df) || !is.data.frame(plant_summaries_df) ||
      nrow(plant_summaries_df) < 1L) {
    return(FALSE)
  }
  if (length(.sf_dr_pc_names(names(plant_with_pcs))) < 1L) {
    return(FALSE)
  }
  identical(
    .sf_dr_plant_fingerprint(plant_summaries_df),
    .sf_dr_plant_fingerprint(plant_with_pcs)
  )
}

#' Prefer PC-augmented plants only while aligned with current summaries.
#' @param plant_summaries_df Current plant summaries (or NULL)
#' @param plant_with_pcs Cached PCA result frame (or NULL)
#' @return data.frame to feed BLUP
.sf_dr_plants_for_blup <- function(plant_summaries_df, plant_with_pcs) {
  if (.sf_dr_pcs_aligned(plant_summaries_df, plant_with_pcs)) {
    return(plant_with_pcs)
  }
  plant_summaries_df
}

.sf_dr_fail <- function(reason, plant_df = NULL, meta_extra = list()) {
  empty_loadings <- data.frame()
  empty_var <- data.frame(
    component = character(0),
    std_dev = numeric(0),
    prop_var = numeric(0),
    cum_prop = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_cluster <- data.frame(
    trait = character(0),
    cluster = integer(0),
    stringsAsFactors = FALSE
  )
  list(
    plant_df_with_pcs = plant_df,
    loadings = empty_loadings,
    variance = empty_var,
    cluster = empty_cluster,
    meta = c(
      list(
        success = FALSE,
        reason = reason,
        trait_cols = character(0),
        n_comp = 0L,
        n_plants = if (is.data.frame(plant_df)) nrow(plant_df) else 0L,
        n_complete = 0L,
        scale = NA
      ),
      meta_extra
    )
  )
}

#' Trait correlation clustering (absolute correlation distance).
#' @param X Complete numeric matrix (plants × traits)
#' @return data.frame(trait, cluster)
.sf_dr_cluster_traits <- function(X) {
  traits <- colnames(X)
  if (length(traits) < 2L) {
    return(data.frame(
      trait = traits,
      cluster = if (length(traits) == 1L) 1L else integer(0),
      stringsAsFactors = FALSE
    ))
  }
  cormat <- stats::cor(X, use = "pairwise.complete.obs")
  cormat[!is.finite(cormat)] <- 0
  diag(cormat) <- 1
  d <- stats::as.dist(1 - abs(cormat))
  hc <- stats::hclust(d, method = "average")
  k <- max(2L, min(4L, length(traits) - 1L))
  cl <- stats::cutree(hc, k = k)
  data.frame(
    trait = traits,
    cluster = as.integer(cl[traits]),
    stringsAsFactors = FALSE
  )
}

#' Plant × traits PCA; attach PC scores to the same plant rows.
#'
#' Never uses genotype-aggregated rows. Design factors (cultivar, treatment,
#' …) are stripped from `trait_cols` even if passed.
#'
#' @param plant_df Plant-level wide data frame (one row per plant)
#' @param trait_cols Phenotype columns; default = all eligible numeric
#' @param n_comp Number of components to keep; NULL = all available
#' @param scale Scale variables before PCA (prcomp scale.)
#' @param require_named_temporal_metrics If TRUE, refuse unlabeled raw phenotype
#'   names (forces GS / trait__metric columns for temporal data).
#' @return list(plant_df_with_pcs, loadings, variance, cluster, meta)
reduceTraitsPCA <- function(plant_df, trait_cols = NULL, n_comp = NULL,
                            scale = TRUE,
                            require_named_temporal_metrics = FALSE) {
  if (is.null(plant_df) || !is.data.frame(plant_df) || nrow(plant_df) < 1L) {
    return(.sf_dr_fail("plant_df must be a non-empty data frame"))
  }

  eligible <- .sf_dr_eligible_traits(plant_df)
  if (is.null(trait_cols)) {
    trait_cols <- eligible
  } else {
    trait_cols <- as.character(trait_cols)
    # Strip design / non-numeric / missing — never PCA on factors
    dropped <- setdiff(trait_cols, eligible)
    trait_cols <- intersect(trait_cols, eligible)
    if (length(dropped) > 0L) {
      # kept in meta below
    }
  }

  if (length(trait_cols) < 2L) {
    return(.sf_dr_fail(
      "Need at least 2 numeric phenotype columns for plant×trait PCA.",
      plant_df = plant_df,
      meta_extra = list(requested_dropped = setdiff(
        if (is.null(trait_cols)) character(0) else trait_cols,
        eligible
      ))
    ))
  }

  if (isTRUE(require_named_temporal_metrics)) {
    ok_nm <- if (exists(".sf_is_named_temporal_metric_col", mode = "function")) {
      vapply(trait_cols, .sf_is_named_temporal_metric_col, logical(1))
    } else {
      grepl("__(auc|peak|growth_rate|last|mean)$", trait_cols) |
        trait_cols %in% c("mean", "auc", "peak", "growth_rate", "last")
    }
    bad <- trait_cols[!ok_nm]
    if (length(bad) > 0L) {
      return(.sf_dr_fail(
        paste0(
          "Named temporal metrics required (trait__auc/peak/growth_rate/last); ",
          "unlabeled raw phenotypes forbidden as silent across-time means: ",
          paste(utils::head(bad, 8L), collapse = ", ")
        ),
        plant_df = plant_df
      ))
    }
  }

  # Drop zero-variance traits (after complete-case subset)
  X_all <- as.matrix(plant_df[, trait_cols, drop = FALSE])
  storage.mode(X_all) <- "double"
  ok <- stats::complete.cases(X_all)
  n_complete <- sum(ok)
  if (n_complete < 3L) {
    return(.sf_dr_fail(
      "Need at least 3 plants with complete trait values for PCA.",
      plant_df = plant_df,
      meta_extra = list(trait_cols = trait_cols, n_complete = n_complete)
    ))
  }

  X <- X_all[ok, , drop = FALSE]
  sds <- apply(X, 2L, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  if (sum(keep) < 2L) {
    return(.sf_dr_fail(
      "Fewer than 2 traits with positive variance after complete-case filter.",
      plant_df = plant_df,
      meta_extra = list(trait_cols = trait_cols, n_complete = n_complete)
    ))
  }
  dropped_zv <- trait_cols[!keep]
  trait_cols <- trait_cols[keep]
  X <- X[, keep, drop = FALSE]

  max_comp <- min(ncol(X), nrow(X) - 1L)
  if (max_comp < 1L) {
    return(.sf_dr_fail(
      "Not enough rank for PCA.",
      plant_df = plant_df,
      meta_extra = list(trait_cols = trait_cols, n_complete = n_complete)
    ))
  }

  if (is.null(n_comp) || is.na(as.integer(n_comp)[1L])) {
    n_comp <- max_comp
  } else {
    n_comp <- as.integer(n_comp)[1L]
    if (n_comp < 1L) {
      n_comp <- 1L
    }
    n_comp <- min(n_comp, max_comp)
  }

  pca <- stats::prcomp(X, center = TRUE, scale. = isTRUE(scale), rank. = n_comp)

  # Full-length scores aligned to plant rows (NA if incomplete)
  score_mat <- matrix(NA_real_, nrow = nrow(plant_df), ncol = n_comp)
  colnames(score_mat) <- paste0("PC", seq_len(n_comp))
  score_mat[ok, ] <- pca$x[, seq_len(n_comp), drop = FALSE]

  # Strip prior PC columns before binding fresh ones
  out_df <- plant_df
  old_pcs <- .sf_dr_pc_names(names(out_df))
  if (length(old_pcs) > 0L) {
    out_df <- out_df[, setdiff(names(out_df), old_pcs), drop = FALSE]
  }
  for (j in seq_len(n_comp)) {
    out_df[[colnames(score_mat)[j]]] <- score_mat[, j]
  }

  rot <- as.matrix(pca$rotation[, seq_len(n_comp), drop = FALSE])
  loadings <- as.data.frame(rot, stringsAsFactors = FALSE)
  loadings$trait <- rownames(rot)
  rownames(loadings) <- NULL
  loadings <- loadings[, c("trait", colnames(score_mat)), drop = FALSE]

  sdev <- as.numeric(pca$sdev[seq_len(n_comp)])
  var_exp <- sdev^2
  total <- sum(pca$sdev^2)
  prop <- if (total > 0) var_exp / total else rep(NA_real_, n_comp)
  variance <- data.frame(
    component = colnames(score_mat),
    std_dev = sdev,
    prop_var = prop,
    cum_prop = cumsum(prop),
    stringsAsFactors = FALSE
  )

  cluster <- .sf_dr_cluster_traits(X)

  list(
    plant_df_with_pcs = out_df,
    loadings = loadings,
    variance = variance,
    cluster = cluster,
    meta = list(
      success = TRUE,
      reason = NA_character_,
      trait_cols = trait_cols,
      dropped_zero_var = dropped_zv,
      n_comp = n_comp,
      n_plants = nrow(plant_df),
      n_complete = n_complete,
      scale = isTRUE(scale)
    )
  )
}
