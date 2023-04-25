

# cleaning -----
rm(list = ls())
gc(reset = TRUE)

# renv -----
renv::activate()
renv::hydrate(prompt = FALSE)
here::i_am('README.md')
set.seed(42)


# project selection -------------------------------------------------------

project <- "project_NO3"

phenospex_file <- Sys.glob(paste0("data/",project, '/*_data.zip'))
unit_file_1 <- Sys.glob(paste0("data/",project, '/*_handmade.csv'))
unit_file_2 <- Sys.glob(paste0("data/",project, '/*_translation.csv'))
group_file <- Sys.glob(paste0("data/",project, '/groups.xlsx'))

remove('times')

# functions -----
`%>%` <- magrittr::`%>%`

sym_diff <- function(a, b) {
  sort(setdiff(union(a, b), intersect(a, b)))
}

tidy_CL_cell <- function(CL_column) {
  #transforms column with objects from multcompView::multcompLetters4(...)
  #into column with tibbles
  CL_column %>%
    purrr::map2(., names(.), ~ {
      tibble::as_tibble(.x$Letters,
                        rownames = 'group') %>%
        dplyr::rename(!!.y := group, tukey_letter = value)
    }) %>%
    dplyr::bind_rows() %>%
    janitor::clean_names() %>%
    tidyr::pivot_longer(
      cols = -tukey_letter,
      names_to = 'tukey_grouping',
      values_to = 'tukey_group'
    ) %>%
    tidyr::drop_na()
}

p_stars <- function(x) {
  #formating of p-values
  dplyr::case_when(x < 0.001 ~ '***',
                   x < 0.01 ~ '**',
                   x < 0.05 ~ '*',
                   x < 0.1 ~ '.',
                   .default = ' ')
}


try_normalize_or_NA <- function(data){
  tryCatch(
    expr = {
      new_val <- bestNormalize::bestNormalize(data,
                                              allow_orderNorm = FALSE,
                                              standardize = FALSE,
                                              out_of_sample = FALSE)$x.t
      return(new_val)
    },
    error = function(e){
      print(e)
      message('bestNormalize returned an error! Deleting the cell')
      return(NA)
    }
  )
}

# import planteye table -----
planteye_table <- readr::read_csv(phenospex_file) %>%
  janitor::clean_names(.)

remove(phenospex_file)

# aggregation with dbscan -----
hours_eps <- 1

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
  dplyr::mutate(dplyr::across(dplyr::contains('_logit'), \(x) dplyr::if_else(is.infinite(x), minus_inf_replacement, x)))

logit_numeric_colnames <- planteye_table %>%
  dplyr::select(dplyr::where(is.numeric)) %>%
  colnames()

# # normalize distribution in each trait_value -----
# planteye_table <- planteye_table %>% # normalize for each trait
#   dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) try_normalize_or_NA(x)))

# import unit data -----
unit_table_1 <- readr::read_csv(unit_file_1)
unit_table_2 <- readr::read_csv(unit_file_2)

unit_table <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = 'V.T.R') %>%
  janitor::clean_names() %>%
  dplyr::rename(repetition = `repeat`)

remove(list = c('unit_file_1', 'unit_file_2',
                'unit_table_1', 'unit_table_2')
       )

# import groups table -----
excel <- readxl::read_xlsx(group_file)

sheetname_vector <- readxl::excel_sheets(group_file) %>%
  unlist()

sheet_list <- sheetname_vector %>%
  purrr::map(\(x) readxl::read_xlsx(group_file, sheet = x))

vector_of_groups <- janitor::make_clean_names(sheetname_vector)
vector_of_groups %>% saveRDS('.cache/vector_of_groups.rds')

groups_table <-
  purrr::map2_dfr(sheetname_vector,
                  sheet_list,
                  \(naming, subtable) subtable %>%
                    dplyr::mutate(source = naming))

groups_table <- groups_table %>%
  tidyr::pivot_longer(contains('Вариант'),
                      names_to = 'group_number',
                      values_to = 'var') %>%
  dplyr::select(!contains('...')) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(group_number = stringi::stri_replace_all_regex(group_number, '[^\\d]', '')) %>%
  dplyr::rename(group_principle = source) %>%
  tidyr::pivot_wider(names_from = group_principle,
                     values_from = group_number) %>%
  janitor::clean_names(.)

remove(list = c('group_file', 'excel',
                'sheet_list', 'sheetname_vector')
       )

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

unit_keys <- unit_table %>%
  dplyr::select(var) %>%
  dplyr::distinct() %>%
  dplyr::pull() %>%
  sort()
group_keys <- groups_table %>%
  dplyr::select(var) %>%
  dplyr::distinct() %>%
  dplyr::pull() %>%
  sort()
grouping_mismatch_keys <- sym_diff(unit_keys, group_keys)

merged_table <- merged_table %>%
  dplyr::inner_join(groups_table, by = 'var')

merged_table <- merged_table %>%
  dplyr::select(-treatment.x) %>%
  dplyr::rename(treatment = treatment.y) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = 'UTC')) %>%
  dplyr::mutate(
    treatment = as.character(treatment),
    repetition = as.character(repetition),
    plant = as.character(plant)
  )

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
  dplyr::select(-total_numeric)  %>%
  tidyr::drop_na(vector_of_groups)

saveRDS(merged_table, file = '.cache/merged_table.rds')

remove(list = ls())
gc(reset = TRUE)


# shiny -------------------------------------------------------------------


shiny::runApp('shiny_experimental', launch.browser = TRUE)
