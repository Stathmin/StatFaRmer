suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(emmeans)
  library(multcompView)
})

cat("Starting CLD debug script...\n")

# Load data
df_path <- here('data','project_soy_2024-05','project_soy_2024-05_merged_table.rds')
if (!file.exists(df_path)) {
  stop(paste('Data not found at', df_path))
}
df <- readRDS(df_path)

# Ensure factor types are consistent
if (!('dbscan_cluster' %in% names(df))) stop('dbscan_cluster missing in data')
if (!('treatment' %in% names(df))) stop('treatment missing in data')
df$dbscan_cluster <- as.factor(df$dbscan_cluster)
df$treatment <- as.factor(df$treatment)

# Filter clusters 1,4,6 as per logs and droplevels
keep_clusters <- intersect(levels(df$dbscan_cluster), c('1','4','6'))
if (length(keep_clusters) == 0) {
  # In case levels are numeric-like
  keep_clusters <- intersect(as.character(sort(unique(df$dbscan_cluster))), c('1','4','6'))
}
df <- df %>% filter(as.character(dbscan_cluster) %in% keep_clusters)
df <- droplevels(df)

outcome <- 'canopy_light_penetration_depth_mm'
if (!(outcome %in% names(df))) stop(paste('Outcome', outcome, 'not in data'))

form <- as.formula(paste(outcome, '~ 1 + (dbscan_cluster + treatment)^2'))
cat('Fitting model:', deparse(form), '\n')
m <- aov(form, data = df)

check_one <- function(f, by) {
  cat('\n=== Factor:', f, 'By:', paste(by, collapse=','), '===\n')
  emm <- if (length(by)) emmeans(m, specs = f, by = by) else emmeans(m, specs = f)
  print(utils::head(summary(emm)))
  pw <- contrast(emm, 'pairwise', adjust='tukey')
  pw_df <- as.data.frame(pw)
  print(utils::head(pw_df))
  na_cnt <- sum(is.na(pw_df$p.value)) + if ('p.adj' %in% names(pw_df)) sum(is.na(pw_df$p.adj)) else 0
  cat('NA p-values count:', na_cnt, '\n')

  cld_tab <- tryCatch(emmeans::cld(emm, adjust='tukey'), error=function(e) NULL)
  if (!is.null(cld_tab)) {
    cols <- c(by, f, '.group')
    cols <- cols[cols %in% names(cld_tab)]
    print(utils::head(cld_tab[, cols, drop=FALSE]))
    blank_cnt <- sum(trimws(cld_tab$.group)=='' | is.na(cld_tab$.group))
    cat('Blank .group count:', blank_cnt, 'of', nrow(cld_tab), '\n')
  } else {
    cat('cld() returned NULL\n')
  }

  if (is.null(cld_tab) || all(trimws(cld_tab$.group)=='')) {
    pvals <- if ('p.adj' %in% names(pw_df)) pw_df$p.adj else pw_df$p.value
    names(pvals) <- pw_df$contrast
    letters <- tryCatch(multcompView::multcompLetters(pvals), error=function(e) NULL)
    cat('multcompLetters result:\n')
    print(letters)
  }

  data_levels <- unique(as.character(df[[f]]))
  cld_levels <- if (!is.null(cld_tab)) unique(as.character(cld_tab[[f]])) else character(0)
  cat('Data levels sample:', paste(utils::head(data_levels), collapse=', '), '\n')
  cat('CLD levels sample:', paste(utils::head(cld_levels), collapse=', '), '\n')
}

check_one('dbscan_cluster', by = c('treatment'))
check_one('treatment', by = c('dbscan_cluster'))

cat('\nDone.\n')

