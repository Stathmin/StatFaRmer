suppressPackageStartupMessages({
	library(here)
	library(emmeans)
})

source(here('shiny','common','anova.R'))
source(here('shiny','common','tukey.R'))
source(here('shiny','common','letters_pvals.R'))

df <- readRDS(here('data','project_soy_2024-05','project_soy_2024-05_merged_table.rds'))
df$dbscan_cluster <- factor(df$dbscan_cluster)
df$treatment <- factor(df$treatment)
keep <- intersect(levels(df$dbscan_cluster), c('1','4','6'))
if (length(keep) > 0) {
	df <- subset(df, as.character(dbscan_cluster) %in% keep)
	df <- droplevels(df)
}

an <- performANOVA(df, 'canopy_light_penetration_depth_mm ~ 1 + (dbscan_cluster + treatment)^2', 'canopy_light_penetration_depth_mm')
stopifnot(an$success)

emm <- emmeans(an$model, specs = 'dbscan_cluster', by = 'treatment')
pw <- as.data.frame(contrast(emm, method = 'pairwise', adjust = 'tukey'))
cat('PW names:', paste(names(pw), collapse=', '), '\n')
print(utils::head(pw))
str(pw)

cat('multcompLetters on pw p.values (named by contrast) -> names/values:\n')
ml <- multcompView::multcompLetters(setNames(pw$p.value, as.character(pw$contrast)))
print(ml)
cat('names(ml$Letters):', paste(names(ml$Letters), collapse=', '), '\n')

cat('buildLettersFromPvals result for dbscan_cluster by treatment:\n')
print(buildLettersFromPvals(pw, 'dbscan_cluster', byNames = c('treatment')))


