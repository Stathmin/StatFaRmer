# cleaning -----
rm(list = ls())
gc(reset = TRUE)

# renv -----
set.seed(42)


# project selection -------------------------------------------------------

project <- "project_NO3" # name of your project subfolder in data folder

hours_eps <- 1 #time clustering distance in hours

phenospex_file <- Sys.glob(paste0("data/",project, '/*_data.zip'))
unit_file_1 <- Sys.glob(paste0("data/",project, '/*_handmade.csv'))
unit_file_2 <- Sys.glob(paste0("data/",project, '/*_translation.csv'))
group_file <- Sys.glob(paste0("data/",project, '/groups.xlsx'))

# functions -----
`%>%` <- magrittr::`%>%`

sym_diff <- function(a, b) {
  sort(setdiff(union(a, b), intersect(a, b)))
}

p_stars <- function(x) {
  #formating of p-values
  dplyr::case_when(x < 0.001 ~ '***',
                   x < 0.01 ~ '**',
                   x < 0.05 ~ '*',
                   x < 0.1 ~ '.',
                   .default = ' ')
}

# import planteye table -----
planteye_table <- readr::read_csv(phenospex_file) %>%
  janitor::clean_names(.)

remove(phenospex_file)

# aggregation with dbscan -----

planteye_table <- planteye_table %>%
  dplyr::mutate(hours_from_start =
                  as.numeric(difftime(timestamp, min(timestamp),
                                      units = 'hours')))

dbscan_cluster <- planteye_table %>%
  dplyr::select(hours_from_start) %>%
  dplyr::pull() %>%
  matrix(ncol = 1) %>%
  dbscan::dbscan(., eps = hours_eps) %>%
  .$cluster %>%
  forcats::as_factor(.) %>%
  tibble::as_tibble_col(column_name = 'dbscan_cluster')

checkmate::assert_true(dbscan_cluster %>%
                         dplyr::n_distinct() > 1)

planteye_table <- planteye_table %>%
  dplyr::bind_cols(dbscan_cluster) %>%
  dplyr::select(-hours_from_start)

remove(list = c('dbscan_cluster', 'hours_eps'))

# remove outlier groups -----

within_groups_not_outliers <- planteye_table %>% #delete over-3-sigmas
  dplyr::group_by(dbscan_cluster) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                              {\(x) (x - min(x))/(max(x) - min(x))})) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                              {\(x) abs(x - mean(x)) / sd(x) <= 3})) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.logical),
                              ~dplyr::if_else(is.na(.x), TRUE, .x))) %>%
  tidyr::pivot_longer(dplyr::where(is.logical),
                      names_to = 'trait',
                      values_to = 'trait_value') %>%
  dplyr::filter(trait_value) %>%
  dplyr::select(-trait_value)


planteye_table <- planteye_table %>%
  tidyr::pivot_longer(dplyr::where(is.numeric),
                      names_to = 'trait',
                      values_to = 'trait_value') %>%
  dplyr::right_join(within_groups_not_outliers) %>%
  tidyr::pivot_wider(names_from = 'trait',
                     values_from = 'trait_value')

remove(within_groups_not_outliers)

# percentage to logit -----
fix_perc_imprecision <- \(x) dplyr::case_when((x >= 0) &
                                                (x <= 1.00) ~ x,
                                              (x < 0) &
                                                (x >= -0.01) ~ 0,
                                              (x > 1) &
                                                (x <= 1.01) ~ 1,
                                              .default = NA)

planteye_table <- planteye_table %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::contains('_percent'),
      \(x) stats::qlogis(fix_perc_imprecision(x))
      )
    ) %>%
  dplyr::rename_with(\(x) stringi::stri_replace_all_fixed(x,
                                                          pattern = '_percent',
                                                          replacement = '_logit'))

# will replace -Inf (from 0 inputs) with the closest negative value
minus_inf_replacement <- planteye_table %>%
  dplyr::select(dplyr::contains('_logit')) %>% dplyr::pull() %>%
  .[!is.na(.) & !is.infinite(.)] %>%
  min() %>%
  floor()

planteye_table <- planteye_table %>%
  dplyr::mutate(dplyr::across(dplyr::contains('_logit'),
                              \(x) dplyr::if_else(is.infinite(x),
                                                  minus_inf_replacement,
                                                  x)))

logit_numeric_colnames <- planteye_table %>%
  dplyr::select(dplyr::where(is.numeric)) %>%
  colnames()

# import unit data -----
unit_table_1 <- readr::read_csv(unit_file_1) %>%
  #dplyr::rename(dplyr::all_of(c(repetition = 'Repeat'))) %>%
  dplyr::mutate_all(as.character)
unit_table_2 <- readr::read_csv(unit_file_2) %>%
  dplyr::mutate_all(as.character)

unit_table <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = 'V.T.R') %>%
  janitor::clean_names()

remove(list = c('unit_file_1', 'unit_file_2',
                'unit_table_1', 'unit_table_2')
)

# import groups table -----
if (length(group_file) > 0) {
  groups_table = readxl::read_xlsx(group_file) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

  remove(list = c('group_file')
         )
}

# table merge -----
dim_initial <- dim(planteye_table)

planteye_keys <- planteye_table %>%
  dplyr::select(unit) %>%
  dplyr::distinct() %>%
  dplyr::pull() %>%
  sort()
unit_keys <- unit_table %>%
  dplyr::select(t_x_y) %>%
  dplyr::distinct() %>%
  dplyr::pull() %>%
  sort()
unit_mismatch_keys <- sym_diff(unit_keys, planteye_keys)

merged_table <- planteye_table %>%
  dplyr::inner_join(unit_table, by = c('unit' = 't_x_y'))

dim_after_translate <- dim(merged_table)

if (exists("groups_table")){
  unit_keys <- unit_table %>%
    dplyr::select(cultivar) %>%
    dplyr::distinct() %>%
    dplyr::pull() %>%
    sort()

  group_keys <- groups_table %>%
    dplyr::select(cultivar) %>%
    dplyr::distinct() %>%
    dplyr::pull() %>%
    sort()

  grouping_mismatch_keys <- sym_diff(unit_keys, group_keys)

  merged_table <- merged_table %>%
    dplyr::inner_join(groups_table, by = 'cultivar')
}

merged_table <- merged_table %>%
  dplyr::select(-treatment.x) %>%
  dplyr::rename(treatment = treatment.y) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = 'UTC')) %>%
  dplyr::mutate(
    treatment = as.character(treatment)
  ) %>%
  janitor::remove_constant()

dim_after_group <- dim(merged_table)

lost_rows <- (dim_initial - dim_after_group)[1]

string_colnames <- merged_table %>%
  dplyr::select(where(\(x) !is.numeric(x) &
                        !lubridate::is.timepoint(x))) %>%
  colnames()

remove(list = c('planteye_table', 'unit_table',
                'groups_table', 'planteye_keys',
                'unit_keys', 'group_keys')
       )

merged_table <- merged_table %>%
  dplyr::relocate(all_of(string_colnames), .before = timestamp) %>%
  dplyr::mutate(total_numeric = rowSums(pick(where(is.numeric)),
                                        na.rm = TRUE)) %>%
  dplyr::filter(total_numeric != 0) %>%
  dplyr::select(-total_numeric)

mean_table <- merged_table %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::where(is.numeric),
      \(x) median(x, na.rm = TRUE)
    )
  ) %>%
  dplyr::ungroup() %>%
  janitor::remove_constant(na.rm = TRUE)

merged_table <- merged_table %>%
  dplyr::select(
    -dplyr::where(is.numeric)
  ) %>%
  dplyr::left_join(mean_table, by = c("unit", "dbscan_cluster")) %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


saveRDS(merged_table, file = 'shiny/merged_table.rds')
merged_table %>%
  dplyr::select(where(is.character)) %>%
  colnames() %>%
  {.[!. %in% c("unit","v_t_r","cultivar")]} %>%
  saveRDS(file = 'shiny/vector_of_groups.rds')

remove(list = ls())
gc(reset = TRUE)


# shiny -------------------------------------------------------------------

#merged_table <- readRDS('shiny/merged_table.rds') #For debug

shiny::runApp('shiny', launch.browser = TRUE)
