# Biological interpretation of Trait-DR PCs + BLUP ranks (and optional raw traits).
# Produces human-readable narratives: axes meaning, best/worst genotypes,
# cross-PC overlaps, and biomass↔NDVI style trade-off groups.
#
# !!! THIS CODE SHOULD BE CHECKED BY A HUMAN PLZ !!!
# Auto-generated PC meanings, GLI-bin polarity, and genotype archetypes are aids —
# not publication-ready conclusions. A human must review before citing.

if (!exists(".sf_hortcontrol_trait_label", mode = "function")) {
  .hc <- file.path("shiny", "common", "hortcontrol_labels.R")
  if (!file.exists(.hc) && requireNamespace("here", quietly = TRUE)) {
    .hc <- here::here("shiny", "common", "hortcontrol_labels.R")
  }
  if (file.exists(.hc)) source(.hc, local = FALSE)
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Pretty trait label for biology text (HortControl-style).
.sf_biology_trait_label <- function(trait_name) {
  .sf_hortcontrol_trait_label(trait_name)
}

#' Map PlantEye / phenotype column names to short biology tags.
.sf_biology_trait_tag <- function(trait_name) {
  nm <- tolower(as.character(trait_name))
  # Strip GS metric suffix (trait__auc → trait)
  base <- sub("__(auc|peak|growth_rate|last|mean)$", "", nm)
  fam <- if (grepl("digital_biomass|biomass|volume", base)) {
    "biomass/size"
  } else if (grepl("leaf_area|3d_leaf|x3d_leaf", base)) {
    "leaf area"
  } else if (grepl("convex_hull|hull_|circumference|width|height|depth|plant_height", base)) {
    "canopy geometry"
  } else if (grepl("^ndvi_.*_percent$", base) || grepl("^ndvi_average$", base)) {
    "NDVI / greenness"
  } else if (grepl("^gli_.*_percent$", base) || grepl("^gli_average$", base)) {
    "GLI / greenness"
  } else if (grepl("saturation", base)) {
    "color saturation"
  } else if (grepl("hue", base)) {
    "hue / color"
  } else if (base %in% c("auc", "peak", "mean", "growth_rate", "last")) {
    "growth-curve summary"
  } else {
    "other phenotype"
  }
  fam
}

#' Whether the high end of a PC means greener tissue (HortControl sense).
#' Anchors on NDVI/GLI *average* loadings; low-index bins are anti-green.
.sf_biology_high_is_greener <- function(vals) {
  nm <- names(vals)
  base <- sub("__(auc|peak|growth_rate|last|mean)$", "", tolower(nm))
  avg <- vals[base %in% c("ndvi_average", "gli_average")]
  if (length(avg) >= 1L && any(is.finite(avg))) {
    return(isTRUE(mean(avg, na.rm = TRUE) > 0))
  }
  # Fallback: high positive bins (upper half of index range) vs low/negative bins
  is_bin <- grepl("percent$", base) & grepl("^(gli|ndvi)_", base)
  if (!any(is_bin)) return(NA)
  iv <- vapply(base[is_bin], .sf_spectral_bin_interval, character(1))
  # Midpoint of interval string "[a, b]"
  mid <- vapply(iv, function(z) {
    if (is.na(z)) return(NA_real_)
    nums <- as.numeric(unlist(regmatches(z, gregexpr("-?[0-9]+(?:\\.[0-9]+)?", z))))
    if (length(nums) < 2L) return(NA_real_)
    mean(nums[1:2])
  }, numeric(1))
  ok <- is.finite(mid)
  if (!any(ok)) return(NA)
  cor(mid[ok], as.numeric(vals[is_bin][ok]), use = "complete.obs") > 0
}

#' Biological axis summary from PCA loadings (signed).
#' @return list(title, blurb, top_pos, top_neg, family)
.sf_biology_axis_from_loadings <- function(loadings, pc, top_n = 5L) {
  empty <- list(
    title = pc, blurb = paste(pc, "not available."),
    top_pos = character(0), top_neg = character(0), family = NA_character_
  )
  if (is.null(loadings) || !is.data.frame(loadings) || !pc %in% names(loadings)) {
    return(empty)
  }
  vals <- setNames(as.numeric(loadings[[pc]]), as.character(loadings$trait))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 1L) return(empty)

  tags <- vapply(names(vals), .sf_biology_trait_tag, character(1))
  fam_score <- tapply(abs(vals), tags, sum, na.rm = TRUE)
  fam_score <- sort(fam_score, decreasing = TRUE)
  family <- names(fam_score)[[1L]]

  # Mean signed loading within dominant family → which end is "more of that"
  fam_traits <- names(vals)[tags == family]
  fam_signed <- vals[fam_traits]
  fam_mean <- mean(fam_signed, na.rm = TRUE)
  high_means_more <- is.finite(fam_mean) && fam_mean > 0
  # Disagreeing signs inside family → contrast axis (not "more green")
  fam_contrast <- length(fam_signed) >= 2L &&
    any(fam_signed > 0.1, na.rm = TRUE) &&
    any(fam_signed < -0.1, na.rm = TRUE)

  high_is_greener <- if (family %in% c("NDVI / greenness", "GLI / greenness")) {
    .sf_biology_high_is_greener(vals)
  } else {
    NA
  }

  ord <- order(vals, decreasing = TRUE)
  top_pos <- names(vals)[ord[seq_len(min(top_n, length(ord)))]]
  top_neg <- names(vals)[rev(ord)[seq_len(min(top_n, length(ord)))]]
  fmt_load <- function(trs) {
    paste(sprintf("%s (%+.2f)", vapply(trs, .sf_biology_trait_label, character(1)), vals[trs]),
          collapse = "; ")
  }

  title <- if (family %in% c("NDVI / greenness", "GLI / greenness")) {
    idx <- if (grepl("^NDVI", family)) "NDVI" else "GLI"
    if (isTRUE(fam_contrast)) {
      if (isTRUE(high_is_greener)) {
        paste0(idx, " distribution: high PC = greener (higher average / high-index bins)")
      } else if (identical(high_is_greener, FALSE)) {
        paste0(idx, " distribution: high PC = more low/negative-index tissue (bins), lower ",
               idx, " average")
      } else {
        paste0(idx, " bin% ↔ average contrast (not a single greenness score)")
      }
    } else if (isTRUE(high_is_greener) || (is.na(high_is_greener) && high_means_more)) {
      paste0(idx, " greenness (+)")
    } else {
      paste0(idx, " greenness (−)")
    }
  } else {
    switch(
      family,
      "biomass/size" = if (high_means_more) "Canopy size / biomass (+)" else "Canopy size / biomass (−)",
      "leaf area" = if (high_means_more) "Leaf area (+)" else "Leaf area (−)",
      "canopy geometry" = if (high_means_more) "Canopy extent / hull (+)" else "Canopy extent / hull (−)",
      "color saturation" = if (high_means_more) "Color saturation (+)" else "Color saturation (−) ↔ greenness mix",
      "hue / color" = "Hue / canopy color axis",
      paste(family, "axis")
    )
  }

  blurb <- sprintf(
    paste0(
      "%s is mainly a **%s** axis. ",
      "Spectral bins are %% of 3D points in an index interval (HortControl); ",
      "averages are mean per-point index. High %s scores load toward: %s. ",
      "Low scores toward: %s."
    ),
    pc, family, pc, fmt_load(top_pos), fmt_load(top_neg)
  )

  list(
    title = title,
    blurb = blurb,
    family = family,
    high_means_more = high_means_more,
    high_is_greener = high_is_greener,
    top_pos = top_pos,
    top_neg = top_neg,
    prop_hint = NULL
  )
}

#' Rank table helper: top/bottom k levels.
.sf_rank_ends <- function(geno_tbl, k = 10L) {
  if (is.null(geno_tbl) || !is.data.frame(geno_tbl) || nrow(geno_tbl) < 1L) {
    return(list(top = geno_tbl, bottom = geno_tbl))
  }
  df <- geno_tbl
  if (!"rank" %in% names(df)) df$rank <- rank(-df$blup, ties.method = "min")
  df <- df[order(df$rank, df$level), , drop = FALSE]
  k <- max(1L, min(as.integer(k), nrow(df)))
  list(
    top = df[seq_len(k), , drop = FALSE],
    bottom = df[seq.int(nrow(df) - k + 1L, nrow(df)), , drop = FALSE],
    n = nrow(df)
  )
}

#' Join ranks across PCs into one wide table.
.sf_join_pc_ranks <- function(blup_by_pc) {
  pcs <- intersect(c("PC1", "PC2", "PC3", "PC4"), names(blup_by_pc))
  out <- NULL
  for (pc in pcs) {
    bl <- blup_by_pc[[pc]]
    if (is.null(bl) || !isTRUE(bl$success) || is.null(bl$genotype_table)) next
    g <- bl$genotype_table[, c("level", "rank", "blup"), drop = FALSE]
    names(g) <- c("genotype", paste0("rank_", tolower(pc)), paste0("blup_", tolower(pc)))
    out <- if (is.null(out)) g else merge(out, g, by = "genotype", all = TRUE)
  }
  out
}

#' Classify genotypes into size↔greenness archetypes using two PC rank columns.
#' @param size_high_is_bigger if FALSE, high BLUP on size_col means *smaller* canopy
#'   (PC with negative biomass/hull loadings) — ends are flipped for labels.
#' @param green_high_is_greener if FALSE, high BLUP/rank-1 on green_col means *less*
#'   green (e.g. PC loaded on low-GLI bins vs gli_average) — ranks are flipped for labels.
#' @param green_label short name for the green axis in trade-off text (NDVI/GLI/…)
.sf_biology_archetypes <- function(rank_df, size_col, green_col, q = 0.25,
                                   size_high_is_bigger = TRUE,
                                   green_high_is_greener = TRUE,
                                   green_label = "greenness") {
  if (is.null(rank_df) || !is.data.frame(rank_df)) return(list())
  if (!size_col %in% names(rank_df) || !green_col %in% names(rank_df)) return(list())
  df <- rank_df[is.finite(rank_df[[size_col]]) & is.finite(rank_df[[green_col]]), , drop = FALSE]
  if (nrow(df) < 8L) return(list())
  n <- nrow(df)
  hi <- max(1L, floor(n * q))
  # lower rank number = higher BLUP on that PC
  sz_hi <- df[[size_col]] <= hi
  sz_lo <- df[[size_col]] > (n - hi)
  if (isFALSE(size_high_is_bigger)) {
    df$.size_hi <- sz_lo
    df$.size_lo <- sz_hi
  } else {
    df$.size_hi <- sz_hi
    df$.size_lo <- sz_lo
  }
  pc_hi <- df[[green_col]] <= hi
  pc_lo <- df[[green_col]] > (n - hi)
  # Map PC ends → biological greenness ends
  if (isFALSE(green_high_is_greener)) {
    df$.green_hi <- pc_lo
    df$.green_lo <- pc_hi
  } else {
    df$.green_hi <- pc_hi
    df$.green_lo <- pc_lo
  }

  pick <- function(mask, label) {
    sub <- df[mask, , drop = FALSE]
    list(
      label = label,
      n = nrow(sub),
      genotypes = as.character(sub$genotype[order(sub[[size_col]], sub[[green_col]])])
    )
  }
  g <- green_label

  list(
    big_and_green = pick(
      df$.size_hi & df$.green_hi,
      sprintf("Large canopy + strong %s (top quartile both)", g)
    ),
    big_weak_green = pick(
      df$.size_hi & df$.green_lo,
      sprintf("Large canopy but weak %s (size↑ / %s↓ trade-off)", g, g)
    ),
    small_and_green = pick(
      df$.size_lo & df$.green_hi,
      sprintf("Small canopy but strong %s", g)
    ),
    small_weak_green = pick(
      df$.size_lo & df$.green_lo,
      sprintf("Small canopy + weak %s (bottom quartile both)", g)
    )
  )
}

#' Optional single-time raw trait genotype means for narrative (not for Tukey).
.sf_biology_raw_trait_ranks <- function(obs_data, traits, time_value = NULL, group = "cultivar") {
  if (is.null(obs_data) || !is.data.frame(obs_data)) return(list())
  traits <- intersect(traits, names(obs_data))
  if (length(traits) < 1L || !group %in% names(obs_data)) return(list())

  # Prefer helpers from trait_followup if loaded
  if (exists(".sf_slice_single_time", mode = "function")) {
    sl <- .sf_slice_single_time(obs_data, time_value = time_value)
    df <- sl$data
    note <- sl$note
  } else {
    df <- obs_data
    note <- "Full table (no time slice helper)"
  }
  if (is.null(df) || nrow(df) < 1L) return(list())

  out <- list(note = note, traits = list())
  for (tr in traits) {
    x <- suppressWarnings(as.numeric(df[[tr]]))
    g <- as.character(df[[group]])
    ok <- is.finite(x) & !is.na(g)
    if (!any(ok)) next
    m <- tapply(x[ok], g[ok], mean, na.rm = TRUE)
    m <- sort(m, decreasing = TRUE)
    out$traits[[tr]] <- list(
      top = names(m)[seq_len(min(10L, length(m)))],
      bottom = names(m)[seq.int(max(1L, length(m) - 9L), length(m))],
      means = m
    )
  }
  out
}

#' Build full biological readout.
#'
#' @param loadings DR loadings
#' @param variance DR variance table
#' @param blup_by_pc named list PC1/PC2/PC3 → fitBlupModels results
#' @param check check genotype id
#' @param obs_data optional observation table for raw biomass/NDVI notes
#' @param top_k ends for lists
#' @param overlap_k "top k on both" threshold
#' @return list(markdown, bullets, axes, ranks, archetypes, overlaps, raw)
buildBiologyReadout <- function(loadings,
                                variance = NULL,
                                blup_by_pc = NULL,
                                check = NULL,
                                obs_data = NULL,
                                top_k = 10L,
                                overlap_k = 25L,
                                time_value = NULL) {
  pcs <- intersect(c("PC1", "PC2", "PC3"), names(blup_by_pc) %||% character(0))
  if (length(pcs) < 1L && is.data.frame(loadings)) {
    pcs <- intersect(c("PC1", "PC2", "PC3"), names(loadings))
  }

  axes <- lapply(pcs, function(pc) {
    ax <- .sf_biology_axis_from_loadings(loadings, pc, top_n = 5L)
    if (is.data.frame(variance) && "component" %in% names(variance)) {
      pv <- variance$prop_var[match(pc, variance$component)]
      if (length(pv) && is.finite(pv[[1]])) {
        ax$prop_hint <- sprintf("%.0f%% of selected-trait variance", 100 * pv[[1]])
        ax$blurb <- paste0(ax$blurb, " (", ax$prop_hint, ").")
      }
    }
    ax
  })
  names(axes) <- pcs

  # Prefer biomass/leaf for size; prefer NDVI/GLI (not saturation) for greenness
  size_pref <- c("biomass/size", "leaf area", "canopy geometry")
  green_pref <- c("NDVI / greenness", "GLI / greenness", "color saturation", "hue / color")
  size_pc <- NULL
  green_pc <- NULL
  for (fam in size_pref) {
    hit <- pcs[vapply(pcs, function(pc) identical(axes[[pc]]$family, fam), logical(1))]
    if (length(hit)) {
      size_pc <- hit[[1L]]
      break
    }
  }
  for (fam in green_pref) {
    hit <- pcs[vapply(pcs, function(pc) identical(axes[[pc]]$family, fam), logical(1))]
    if (length(hit)) {
      green_pc <- hit[[1L]]
      break
    }
  }
  if (is.null(size_pc) && "PC3" %in% pcs) size_pc <- "PC3"
  if (is.null(green_pc) && "PC2" %in% pcs) green_pc <- "PC2"

  rank_df <- .sf_join_pc_ranks(blup_by_pc)
  ends <- list()
  for (pc in pcs) {
    bl <- blup_by_pc[[pc]]
    if (isTRUE(bl$success)) ends[[pc]] <- .sf_rank_ends(bl$genotype_table, k = top_k)
  }

  overlaps <- list()
  if (!is.null(rank_df) && nrow(rank_df) > 0L) {
    pairs <- list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))
    for (pr in pairs) {
      c1 <- paste0("rank_", tolower(pr[[1]]))
      c2 <- paste0("rank_", tolower(pr[[2]]))
      if (!all(c(c1, c2) %in% names(rank_df))) next
      both <- rank_df[
        is.finite(rank_df[[c1]]) & is.finite(rank_df[[c2]]) &
          rank_df[[c1]] <= overlap_k & rank_df[[c2]] <= overlap_k,
        , drop = FALSE
      ]
      both <- both[order(both[[c1]], both[[c2]]), , drop = FALSE]
      overlaps[[paste(pr, collapse = "_")]] <- both
    }
  }

  archetypes <- list()
  if (!is.null(size_pc) && !is.null(green_pc) && !is.null(rank_df)) {
    g_ax <- axes[[green_pc]]
    s_ax <- axes[[size_pc]]
    high_green <- if (!is.null(g_ax$high_is_greener) && is.logical(g_ax$high_is_greener) &&
                        length(g_ax$high_is_greener) == 1L && !is.na(g_ax$high_is_greener)) {
      isTRUE(g_ax$high_is_greener)
    } else {
      TRUE
    }
    # Canopy/biomass (−) title ⇒ high PC score = smaller plants
    high_bigger <- if (!is.null(s_ax$high_means_more) && is.logical(s_ax$high_means_more) &&
                         length(s_ax$high_means_more) == 1L && !is.na(s_ax$high_means_more)) {
      isTRUE(s_ax$high_means_more)
    } else {
      TRUE
    }
    g_lab <- if (identical(g_ax$family, "NDVI / greenness")) {
      "NDVI"
    } else if (identical(g_ax$family, "GLI / greenness")) {
      "GLI greenness"
    } else {
      "greenness"
    }
    archetypes <- .sf_biology_archetypes(
      rank_df,
      size_col = paste0("rank_", tolower(size_pc)),
      green_col = paste0("rank_", tolower(green_pc)),
      q = 0.25,
      size_high_is_bigger = high_bigger,
      green_high_is_greener = high_green,
      green_label = g_lab
    )
    attr(archetypes, "size_pc") <- size_pc
    attr(archetypes, "green_pc") <- green_pc
    attr(archetypes, "green_high_is_greener") <- high_green
  }

  raw <- .sf_biology_raw_trait_ranks(
    obs_data,
    traits = c("digital_biomass_mm3", "ndvi_average", "plant_height_averaged_mm"),
    time_value = time_value
  )

  # Raw biomass × NDVI archetypes (clearer than PC contrast axes when both exist)
  raw_archetypes <- list()
  if (all(c("digital_biomass_mm3", "ndvi_average") %in% names(raw$traits %||% list()))) {
    bm <- raw$traits$digital_biomass_mm3$means
    nv <- raw$traits$ndvi_average$means
    common <- intersect(names(bm), names(nv))
    if (length(common) >= 8L) {
      rdf <- data.frame(
        genotype = common,
        rank_size = rank(-bm[common], ties.method = "min"),
        rank_green = rank(-nv[common], ties.method = "min"),
        stringsAsFactors = FALSE
      )
      raw_archetypes <- .sf_biology_archetypes(
        rdf, size_col = "rank_size", green_col = "rank_green", q = 0.25
      )
      # Relabel for raw traits
      if (length(raw_archetypes)) {
        raw_archetypes$big_and_green$label <- "High biomass + high NDVI (raw, top quartile both)"
        raw_archetypes$big_weak_green$label <- "High biomass but low NDVI (raw trade-off)"
        raw_archetypes$small_and_green$label <- "Low biomass but high NDVI (raw)"
        raw_archetypes$small_weak_green$label <- "Low biomass + low NDVI (raw)"
        attr(raw_archetypes, "size_pc") <- "digital_biomass_mm3"
        attr(raw_archetypes, "green_pc") <- "ndvi_average"
      }
    }
  }

  bullets <- .sf_biology_bullets(
    axes = axes,
    ends = ends,
    overlaps = overlaps,
    archetypes = archetypes,
    raw_archetypes = raw_archetypes,
    check = check,
    overlap_k = overlap_k,
    size_pc = size_pc,
    green_pc = green_pc,
    raw = raw
  )
  markdown <- .sf_biology_markdown(
    bullets = bullets,
    axes = axes,
    ends = ends,
    overlaps = overlaps,
    archetypes = archetypes,
    raw_archetypes = raw_archetypes,
    check = check,
    overlap_k = overlap_k,
    size_pc = size_pc,
    green_pc = green_pc,
    raw = raw
  )

  list(
    markdown = markdown,
    bullets = bullets,
    axes = axes,
    ends = ends,
    overlaps = overlaps,
    archetypes = archetypes,
    raw_archetypes = raw_archetypes,
    rank_df = rank_df,
    raw = raw,
    size_pc = size_pc,
    green_pc = green_pc,
    check = check
  )
}

.sf_fmt_geno_list <- function(x, max_show = 8L) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) < 1L) return("(none)")
  if (length(x) > max_show) {
    paste0(paste(x[seq_len(max_show)], collapse = ", "), sprintf(" … (+%d)", length(x) - max_show))
  } else {
    paste(x, collapse = ", ")
  }
}

.sf_biology_bullets <- function(axes, ends, overlaps, archetypes, check,
                                overlap_k, size_pc, green_pc, raw,
                                raw_archetypes = list()) {
  out <- character(0)
  for (pc in names(axes)) {
    ax <- axes[[pc]]
    bit <- sprintf("%s ≈ %s", pc, ax$title)
    if (!is.null(ax$prop_hint)) bit <- paste0(bit, " (", ax$prop_hint, ")")
    out <- c(out, bit)
  }
  for (pc in names(ends)) {
    e <- ends[[pc]]
    out <- c(out, sprintf(
      "%s best: %s | worst: %s",
      pc,
      .sf_fmt_geno_list(e$top$level, 5L),
      .sf_fmt_geno_list(rev(e$bottom$level), 5L)
    ))
  }
  key <- paste(c(size_pc, green_pc)[!vapply(c(size_pc, green_pc), is.null, logical(1))], collapse = "_")
  # normalize key order PC2_PC3
  if (!is.null(size_pc) && !is.null(green_pc)) {
    okey <- paste(sort(c(size_pc, green_pc)), collapse = "_")
    # overlaps keys are PC1_PC2 style in declaration order
    for (nm in names(overlaps)) {
      both <- overlaps[[nm]]
      if (is.null(both) || nrow(both) < 1L) {
        out <- c(out, sprintf("Top-%d on both %s: none", overlap_k, gsub("_", " & ", nm)))
      } else {
        out <- c(out, sprintf(
          "Top-%d on both %s (%d gens): %s",
          overlap_k, gsub("_", " & ", nm), nrow(both),
          .sf_fmt_geno_list(both$genotype, 8L)
        ))
      }
    }
  } else {
    for (nm in names(overlaps)) {
      both <- overlaps[[nm]]
      if (!is.null(both) && nrow(both) > 0L) {
        out <- c(out, sprintf(
          "Top-%d on both %s (%d): %s",
          overlap_k, gsub("_", " & ", nm), nrow(both),
          .sf_fmt_geno_list(both$genotype, 8L)
        ))
      }
    }
  }

  if (length(archetypes)) {
    sp <- attr(archetypes, "size_pc") %||% size_pc %||% "size-PC"
    gp <- attr(archetypes, "green_pc") %||% green_pc %||% "green-PC"
    out <- c(out, sprintf("PC trade-off groups (%s × %s, quartiles):", sp, gp))
    for (ar in archetypes) {
      if (is.null(ar) || ar$n < 1L) next
      out <- c(out, sprintf("• %s (n=%d): %s", ar$label, ar$n, .sf_fmt_geno_list(ar$genotypes, 6L)))
    }
  }
  if (length(raw_archetypes)) {
    out <- c(out, "Raw single-time biomass × NDVI groups (quartiles):")
    for (ar in raw_archetypes) {
      if (is.null(ar) || ar$n < 1L) next
      out <- c(out, sprintf("• %s (n=%d): %s", ar$label, ar$n, .sf_fmt_geno_list(ar$genotypes, 6L)))
    }
  }

  if (is.list(raw) && length(raw$traits)) {
    out <- c(out, paste0("Single-time raw means (", raw$note %||% "slice", "):"))
    for (tr in names(raw$traits)) {
      rt <- raw$traits[[tr]]
      out <- c(out, sprintf(
        "• %s top: %s | bottom: %s",
        tr, .sf_fmt_geno_list(rt$top, 5L), .sf_fmt_geno_list(rev(rt$bottom), 5L)
      ))
    }
  }

  if (length(check) && any(nzchar(as.character(check)))) {
    out <- c(out, sprintf("Dunnett check genotype(s): %s", paste(check, collapse = ", ")))
  }
  out
}

.sf_biology_markdown <- function(bullets, axes, ends, overlaps, archetypes,
                                 check, overlap_k, size_pc, green_pc, raw,
                                 raw_archetypes = list()) {
  lines <- c(
    "# Biological readout (StatFaRmer)",
    "",
    "> **!!! THIS CODE SHOULD BE CHECKED BY A HUMAN PLZ !!!**  ",
    "> Auto-labeled PC axes, GLI/NDVI bin polarity, and genotype lists are draft aids — not final biology.",
    "",
    "Interpretation of multi-trait PCA axes and BLUP genotype ranks. ",
    "Ranks: **1 = highest BLUP** on that PC. Single-time raw means are descriptive only ",
    "(not across-time averages; not Tukey letters on 325 genotypes).",
    "",
    "## HortControl spectral parameters (manual §4.3)",
    "",
    "- **GLI / NDVI average:** mean per-point index over all 3D points in the unit.",
    "- **GLI / NDVI bin (%):** percentage of 3D points whose index falls in a named interval ",
    "  (e.g. `GLI bin [0, 0.15]`, `GLI bin [-0.2, 0]`, `GLI bin [-1.0, -0.2]`).",
    "- Low-index bins are **not** “more green” — they are the share of weakly/non-green tissue; ",
    "  greener canopies tend to have higher averages and/or more mass in high-index bins.",
    "",
    "## What the PCs mean biologically",
    ""
  )
  for (pc in names(axes)) {
    ax <- axes[[pc]]
    lines <- c(lines, sprintf("### %s — %s", pc, ax$title), "", ax$blurb, "")
  }

  lines <- c(lines, "## Best / worst genotypes by PC (BLUP)", "")
  for (pc in names(ends)) {
    e <- ends[[pc]]
    lines <- c(
      lines,
      sprintf("### %s", pc),
      "",
      sprintf("- **Top %d:** %s", nrow(e$top), .sf_fmt_geno_list(e$top$level, 15L)),
      sprintf("- **Bottom %d:** %s", nrow(e$bottom), .sf_fmt_geno_list(e$bottom$level, 15L)),
      ""
    )
  }

  lines <- c(lines, sprintf("## Cross-PC overlaps (top %d on both)", overlap_k), "")
  for (nm in names(overlaps)) {
    both <- overlaps[[nm]]
    lab <- gsub("_", " & ", nm)
    if (is.null(both) || nrow(both) < 1L) {
      lines <- c(lines, sprintf("- **%s:** none", lab), "")
    } else {
      lines <- c(
        lines,
        sprintf("- **%s** — %d genotypes: %s", lab, nrow(both), .sf_fmt_geno_list(both$genotype, 20L)),
        ""
      )
    }
  }

  if (length(archetypes)) {
    sp <- attr(archetypes, "size_pc") %||% size_pc
    gp <- attr(archetypes, "green_pc") %||% green_pc
    lines <- c(
      lines,
      sprintf("## PC size ↔ greenness groups (%s × %s, quartiles)", sp, gp),
      ""
    )
    for (ar in archetypes) {
      if (is.null(ar) || ar$n < 1L) {
        lines <- c(lines, sprintf("- **%s:** none", ar$label %||% "?"), "")
        next
      }
      lines <- c(
        lines,
        sprintf("- **%s** (n=%d): %s", ar$label, ar$n, .sf_fmt_geno_list(ar$genotypes, 15L)),
        ""
      )
    }
  }

  if (length(raw_archetypes)) {
    lines <- c(lines, "## Raw biomass × NDVI groups (single time, quartiles)", "")
    for (ar in raw_archetypes) {
      if (is.null(ar) || ar$n < 1L) next
      lines <- c(
        lines,
        sprintf("- **%s** (n=%d): %s", ar$label, ar$n, .sf_fmt_geno_list(ar$genotypes, 15L)),
        ""
      )
    }
  }

  if (is.list(raw) && length(raw$traits)) {
    lines <- c(lines, "## Single-time raw trait leaders", "", raw$note %||% "", "")
    for (tr in names(raw$traits)) {
      rt <- raw$traits[[tr]]
      lines <- c(
        lines,
        sprintf("### `%s`", tr),
        sprintf("- Top: %s", .sf_fmt_geno_list(rt$top, 10L)),
        sprintf("- Bottom: %s", .sf_fmt_geno_list(rev(rt$bottom), 10L)),
        ""
      )
    }
  }

  if (length(check) && any(nzchar(as.character(check)))) {
    lines <- c(lines, "## Check", "", paste(check, collapse = ", "), "")
  }

  lines <- c(lines, "## Short bullets", "", paste0("- ", bullets))
  paste(lines, collapse = "\n")
}

#' Short coach-facing biology paragraph (truncated).
.sf_biology_coach_blurb <- function(readout, max_chars = 900L) {
  if (is.null(readout) || !is.list(readout)) return(NULL)
  b <- readout$bullets
  if (length(b) < 1L) return(NULL)
  txt <- paste(b, collapse = " · ")
  if (nchar(txt) > max_chars) {
    txt <- paste0(substr(txt, 1L, max_chars - 1L), "…")
  }
  txt
}
