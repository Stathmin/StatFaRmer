suppressPackageStartupMessages({
	library(here)
	library(dplyr)
	library(emmeans)
})

cat("Starting CLD soy letters assertion...\n")

source(here('shiny','common','anova.R'))
source(here('shiny','common','tukey.R'))
source(here('shiny','common','letters_utils.R'))
source(here('shiny','common','letters_emmeans.R'))
source(here('shiny','common','letters_pvals.R'))
source(here('shiny','common','letters_join.R'))
source(here('shiny','common','letters_normalize.R'))

df_path <- here('data','project_soy_2024-05','project_soy_2024-05_merged_table.rds')
stopifnot(file.exists(df_path))
df <- readRDS(df_path)

stopifnot('dbscan_cluster' %in% names(df))
stopifnot('treatment' %in% names(df))

df$dbscan_cluster <- factor(df$dbscan_cluster)
df$treatment <- factor(df$treatment)

# Use a variable from logs
outcome <- 'canopy_light_penetration_depth_mm'
stopifnot(outcome %in% names(df))

# Keep a few clusters to match the app case
keep_clusters <- intersect(levels(df$dbscan_cluster), c('1','4','6'))
if (length(keep_clusters) > 0) {
	df <- df %>% filter(as.character(dbscan_cluster) %in% keep_clusters) %>% droplevels()
}

form_text <- paste(outcome, '~ 1 + (dbscan_cluster + treatment)^2')
an <- performANOVA(df, form_text, outcome)
stopifnot(an$success)

tk <- performTukey(an$model, c('dbscan_cluster','treatment'))
stopifnot(tk$success)

letters_list <- lapply(tk$results, function(x) x$letters)
stopifnot(length(letters_list) == 2)

# Assert letters exist and are not all blank/NA
for (i in seq_along(letters_list)) {
	ldf <- letters_list[[i]]
	stopifnot(is.data.frame(ldf))
	stopifnot('group' %in% names(ldf))
	stopifnot('letter' %in% names(ldf))
	all_blank <- all(is.na(ldf$letter) | trimws(ldf$letter) == '')
	if (all_blank) {
		print(ldf)
		stop('Letters are all blank/NA for factor ', names(letters_list)[i])
	}
}

# Exact expectations based on tests/debug_cld.R outputs
# 1) dbscan_cluster within each treatment: 1=a, 4=b, 6=c for var1 and var2
dbscan_letters <- tk$results[["dbscan_cluster"]]$letters
stopifnot(all(c('group','letter','treatment') %in% names(dbscan_letters)))

expected_dbscan <- data.frame(
	group = c('1','4','6','1','4','6'),
	letter = c('a','b','c','a','b','c'),
	treatment = c('var1','var1','var1','var2','var2','var2'),
	stringsAsFactors = FALSE
)
# Normalize types
dbscan_letters$group <- as.character(dbscan_letters$group)
dbscan_letters$treatment <- as.character(dbscan_letters$treatment)
dbscan_letters$letter <- as.character(dbscan_letters$letter)

dbscan_letters <- dbscan_letters[, c('group','letter','treatment')]
dbscan_letters <- dbscan_letters[order(dbscan_letters$treatment, dbscan_letters$group), ]
rownames(dbscan_letters) <- NULL

expected_dbscan <- expected_dbscan[order(expected_dbscan$treatment, expected_dbscan$group), ]
rownames(expected_dbscan) <- NULL

if (!identical(dbscan_letters, expected_dbscan)) {
	print(dbscan_letters)
	print(expected_dbscan)
	stop('dbscan_cluster letters do not match expected per-treatment A/B/C pattern')
}

# 2) treatment within each dbscan_cluster:
# cluster 1: var1=a, var2=a; cluster 4: var1=a, var2=b; cluster 6: var1=a, var2=b
treatment_letters <- tk$results[["treatment"]]$letters
stopifnot(all(c('group','letter','dbscan_cluster') %in% names(treatment_letters)))

expected_treatment <- data.frame(
	group = c('var1','var2','var1','var2','var1','var2'),
	letter = c('a','a','a','b','a','b'),
	dbscan_cluster = c('1','1','4','4','6','6'),
	stringsAsFactors = FALSE
)

treatment_letters$group <- as.character(treatment_letters$group)
treatment_letters$dbscan_cluster <- as.character(treatment_letters$dbscan_cluster)
treatment_letters$letter <- as.character(treatment_letters$letter)

treatment_letters <- treatment_letters[, c('group','letter','dbscan_cluster')]
treatment_letters <- treatment_letters[order(treatment_letters$dbscan_cluster, treatment_letters$group), ]
rownames(treatment_letters) <- NULL

expected_treatment <- expected_treatment[order(expected_treatment$dbscan_cluster, expected_treatment$group), ]
rownames(expected_treatment) <- NULL

if (!identical(treatment_letters, expected_treatment)) {
	print(treatment_letters)
	print(expected_treatment)
	stop('treatment letters do not match expected per-cluster pattern')
}

cat('CLD soy letters assertion passed.\n')


