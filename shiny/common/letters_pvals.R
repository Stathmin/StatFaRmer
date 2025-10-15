buildLettersFromPvals <- function(pw_df, factorName, byNames = NULL) {
  # Fallback CLD per stratum using multcompView::multcompLetters on pw_df
  if (is.null(pw_df) || nrow(pw_df) == 0 || !('contrast' %in% names(pw_df))) {
    return(data.frame(group = character(0), letter = character(0), stringsAsFactors = FALSE))
  }
  # Determine strata
  strata_cols <- if (!is.null(byNames) && length(byNames) > 0) intersect(byNames, names(pw_df)) else character(0)
  if (length(strata_cols) == 0) strata_cols <- character(0)
  # Split by strata (or a single group if none)
  splits <- if (length(strata_cols) > 0) split(pw_df, pw_df[strata_cols], drop = TRUE) else list(`__all__` = pw_df)
  out_list <- list()
  for (nm in names(splits)) {
    sub <- splits[[nm]]
    # Extract and coerce p-values, allowing inputs like "<.0001"
    pvals_raw <- if ('p.adj' %in% names(sub)) sub$p.adj else if ('p.value' %in% names(sub)) sub$p.value else NULL
    if (is.null(pvals_raw)) next
    if (!is.numeric(pvals_raw)) {
      pchar <- as.character(pvals_raw)
      pchar <- sub('^<\\n*\\t*\\r*\\s*', '', pchar)
      pchar <- gsub('[^0-9eE\\.+-]', '', pchar)
      pvals_raw <- suppressWarnings(as.numeric(pchar))
    }
    # Parse contrasts into clean tokens (lhs, rhs)
    contr <- as.character(sub$contrast)
    lr <- strsplit(contr, ' - ', fixed = FALSE)
    lhs <- trimws(vapply(lr, function(p) if (length(p) >= 1) as.character(p[[1]]) else NA_character_, character(1)))
    rhs <- trimws(vapply(lr, function(p) if (length(p) >= 2) as.character(p[[2]]) else NA_character_, character(1)))
    keep <- !(is.na(lhs) | is.na(rhs) | is.na(pvals_raw))
    if (!any(keep)) next
    lhs <- lhs[keep]
    rhs <- rhs[keep]
    pvals <- pvals_raw[keep]
    # Canonicalize pair names so each unordered pair appears once
    # Determine an ordering for tokens (numeric-aware)
    all_tokens <- unique(c(lhs, rhs))
    is_num <- suppressWarnings(!any(is.na(as.numeric(all_tokens))))
    token_order <- if (is_num) order(as.numeric(all_tokens)) else order(all_tokens)
    ordered_levels <- all_tokens[token_order]
    # Build a map from unordered pair -> p-value (prefer minimum if duplicates occur)
    canon_pair <- function(a, b) {
      if (is_num) {
        a_n <- as.numeric(a); b_n <- as.numeric(b)
        if (is.na(a_n) || is.na(b_n)) return(paste(sort(c(a, b)), collapse = '-'))
        if (a_n <= b_n) paste(a, b, sep = '-') else paste(b, a, sep = '-')
      } else {
        if (a <= b) paste(a, b, sep = '-') else paste(b, a, sep = '-')
      }
    }
    pair_names <- mapply(canon_pair, lhs, rhs, USE.NAMES = FALSE)
    # Aggregate duplicates by taking the minimum adjusted p-value (most conservative for grouping)
    agg <- tapply(pvals, pair_names, function(x) min(x, na.rm = TRUE))
    agg <- agg[!is.na(agg)]
    if (length(agg) == 0) next
    # Build a symmetric logical decision matrix using exact level order
    # TRUE = significantly different (per multcompView::multcompLetters logical matrix input)
    alpha <- 0.05
    nlev <- length(ordered_levels)
    M <- matrix(TRUE, nrow = nlev, ncol = nlev, dimnames = list(ordered_levels, ordered_levels))
    if (nlev >= 2) {
      for (i in seq_len(nlev - 1)) {
        for (j in seq.int(i + 1, nlev)) {
          a <- ordered_levels[i]; b <- ordered_levels[j]
          key <- canon_pair(a, b)
          pv <- if (!is.null(agg[[key]]) && is.finite(agg[[key]])) agg[[key]] else NA_real_
          sig <- !is.na(pv) && pv < alpha
          M[a, b] <- sig
          M[b, a] <- sig
        }
      }
    }
    letters_result <- tryCatch(multcompView::multcompLetters(M), error = function(e) NULL)
    if (is.null(letters_result)) next
    # Map letters to exact level order determined above
    letter_map <- as.character(letters_result$Letters)
    names(letter_map) <- trimws(names(letters_result$Letters))
    groups <- ordered_levels
    letters_vec <- as.character(letter_map[groups])
    sub_out <- data.frame(group = groups, letter = letters_vec, stringsAsFactors = FALSE)
    # Add strata columns back
    if (length(strata_cols) > 0 && nm != '__all__') {
      # Parse split name into values when split() formats names as key=value,...
      vals <- lapply(strata_cols, function(sc) unique(as.character(sub[[sc]]))[1])
      for (i in seq_along(strata_cols)) sub_out[[strata_cols[i]]] <- vals[[i]]
    }
    out_list[[length(out_list) + 1]] <- sub_out
  }
  if (length(out_list) == 0) return(data.frame(group = character(0), letter = character(0), stringsAsFactors = FALSE))
  do.call(rbind, out_list)
}


