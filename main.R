# cleaning -----
rm(list = ls())
gc(reset = TRUE)

# renv -----
set.seed(42)


# project selection -----

project <- "project_NO3" # name of your project subfolder in data folder

hours_eps <- 1 #DBSCAN time clustering distance in hours
use_IQR <- FALSE #FALSE -> 3 sigmas; TRUE -> 3 IQR

phenospex_file <- Sys.glob(paste0("data/", project, '/*_data.zip'))
unit_file_1 <- Sys.glob(paste0("data/", project, '/*_handmade.csv'))
unit_file_2 <- Sys.glob(paste0("data/", project, '/*_translation.csv'))
group_file <- Sys.glob(paste0("data/", project, '/groups.xlsx'))

# functions -----
`%>%` <- magrittr::`%>%`

sym_diff <- function(a, b) {
  sort(setdiff(union(a, b), intersect(a, b)))
}

p_stars <- function(x) {
  #formating of p-values
  dplyr::case_when(x < 0.001 ~ '***', x < 0.01 ~ '**', x < 0.05 ~ '*', x < 0.1 ~ '.', .default = ' ')
}

# import planteye table -----
planteye_table <- readr::read_csv(phenospex_file) %>%
  janitor::clean_names(.)

remove(phenospex_file)

# aggregation with dbscan -----

planteye_table <- planteye_table %>%
  dplyr::mutate(hours_from_start =
                  as.numeric(difftime(timestamp, min(timestamp), units = 'hours')))

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

# removal of rows with all observations at zero -----
planteye_table <- planteye_table %>%
  dplyr::filter(rowSums(dplyr::select(., where(is.numeric)) == 0, na.rm = TRUE) <
                  ncol(dplyr::select(., where(is.numeric))))

# percentage to logit transformation -----
fix_perc_imprecision <- \(x) dplyr::case_when((x >= 0) &
                                                (x <= 1.00) ~ x,
                                              (x < 0) &
                                                (x >= -0.01) ~ 0,
                                              (x > 1) &
                                                (x <= 1.01) ~ 1,
                                              .default = NA)

planteye_table <- planteye_table %>%
  dplyr::mutate(dplyr::across(
    dplyr::contains('_percent'),
    \(x) stats::qlogis(fix_perc_imprecision(x))
  )) %>%
  dplyr::rename_with(\(x) stringi::stri_replace_all_fixed(x, pattern = '_percent', replacement = '_logit'))

# replace -Inf after logit (from 0 inputs) with the closest negative value -----
minus_inf_replacement <- planteye_table %>%
  dplyr::select(dplyr::contains('_logit')) %>% dplyr::pull() %>%
  .[!is.na(.) & !is.infinite(.)] %>%
  min() %>%
  floor()

planteye_table <- planteye_table %>%
  dplyr::mutate(dplyr::across(
    dplyr::contains('_logit'),
    \(x) dplyr::if_else(is.infinite(x), minus_inf_replacement, x)
  ))

# remove outlier groups -----

within_groups_not_outliers <- planteye_table %>%
  dplyr::group_by(dbscan_cluster) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), {
    \(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) {
    if (use_IQR) {
      q1 <- quantile(x, 0.25, na.rm = TRUE)
      q3 <- quantile(x, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      return(x >= (q1 - 1.5 * iqr) &
               x <= (q3 + 1.5 * iqr))
    } else {
      return(abs(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) <= 3)
    }
  })) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.logical), ~ dplyr::if_else(is.na(.x), TRUE, .x))) %>%
  tidyr::pivot_longer(dplyr::where(is.logical),
                      names_to = 'trait',
                      values_to = 'trait_value') %>%
  dplyr::filter(trait_value) %>%
  dplyr::select(-trait_value)


normal_table <- planteye_table %>%
  tidyr::pivot_longer(dplyr::where(is.numeric),
                      names_to = 'trait',
                      values_to = 'trait_value') %>%
  dplyr::right_join(within_groups_not_outliers) %>%
  tidyr::pivot_wider(names_from = 'trait', values_from = 'trait_value')

outlier_table <- planteye_table %>%
  tidyr::pivot_longer(dplyr::where(is.numeric),
                      names_to = 'trait',
                      values_to = 'trait_value') %>%
  dplyr::anti_join(within_groups_not_outliers) %>%
  tidyr::pivot_wider(names_from = 'trait', values_from = 'trait_value')

remove(within_groups_not_outliers)

# import unit data -----
unit_table_1 <- readr::read_csv(unit_file_1) %>%
  dplyr::mutate_all(as.character)
unit_table_2 <- readr::read_csv(unit_file_2) %>%
  dplyr::mutate_all(as.character)

unit_table <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = 'V.T.R') %>%
  janitor::clean_names()

remove(list = c('unit_file_1', 'unit_file_2', 'unit_table_1', 'unit_table_2'))

# import groups table -----
if (length(group_file) > 0) {
  groups_table = readxl::read_xlsx(group_file) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

  remove(list = c('group_file'))
}

# table merge -----

all_measurments <- list(normal_table, outlier_table)
remove(list = c('planteye_table', 'normal_table', 'outlier_table'))

all_measurments <- all_measurments %>% purrr::map( ~ {
  initial_df <- .x
  dim_initial <- dim(initial_df)

  planteye_keys <- initial_df %>%
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

  merged_table <- initial_df %>%
    dplyr::inner_join(unit_table, by = c('unit' = 't_x_y'))

  dim_after_translate <- dim(merged_table)

  if (exists("groups_table")) {
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
    dplyr::mutate(treatment = as.character(treatment)) %>%
    janitor::remove_constant()

  dim_after_group <- dim(merged_table)

  lost_rows <- (dim_initial - dim_after_group)[1]

  merged_table

})

merged_table <- all_measurments[[1]]
outlier_table <- all_measurments[[2]]

remove(list = c(
  'unit_table',
  'groups_table',
  'planteye_keys',
  'unit_keys',
  'group_keys'
))

# table reordering -----
string_colnames <- merged_table %>%
dplyr::select(where(\(x) ! is.numeric(x) &
                      !lubridate::is.timepoint(x))) %>%
  colnames()
merged_table <- merged_table %>%
  dplyr::relocate(all_of(string_colnames), .before = timestamp)

# medians for unit technical repetitions within time clusters -----
mean_table <- merged_table %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  janitor::remove_constant(na.rm = TRUE)

merged_table <- merged_table %>%
  dplyr::select(-dplyr::where(is.numeric)) %>%
  dplyr::left_join(mean_table, by = c("unit", "dbscan_cluster")) %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# medians for outliers within time clusters -----
mean_table <- outlier_table %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  janitor::remove_constant(na.rm = TRUE)

outlier_table <- outlier_table %>%
  dplyr::select(-dplyr::where(is.numeric)) %>%
  dplyr::left_join(mean_table, by = c("unit", "dbscan_cluster")) %>%
  dplyr::group_by(unit, dbscan_cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


# exports -----

bind_rows(merged_table %>% mutate(outlier = FALSE),
          outlier_table %>% mutate(outlier = TRUE)) %>%
  mutate()
  saveRDS(file = 'shiny/merged_table.rds')

merged_table %>%
  dplyr::select(where(is.character), -outlier) %>%
  colnames() %>%
  {
    .[!. %in% c("unit", "v_t_r", "cultivar")]
  } %>%
  saveRDS(file = 'shiny/vector_of_groups.rds')

remove(list = ls())
gc(reset = TRUE)


# shiny -----

merged_table <- readRDS('shiny/merged_table.rds') #For debug
vector_of_groups <- readRDS('shiny/vector_of_groups.rds') #For debug

shiny::runApp('shiny', launch.browser = TRUE)
