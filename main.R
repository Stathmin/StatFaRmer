# init --------------------------------------------------------------------
rm(list = ls())
gc()

renv::activate()
renv::hydrate(prompt = FALSE)
here::i_am('README.md')

set.seed(42)
hours_eps <- 0.5
`%>%` <- magrittr::`%>%`

IGNORE_REPORTS <- FALSE

# functions ---------------------------------------------------------------

tidy_CL_cell <- function(CL_column) {
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
  dplyr::case_when(x < 0.001 ~ '***',
                   x < 0.01 ~ '**',
                   x < 0.05 ~ '*',
                   x < 0.1 ~ '.',
                   .default = ' ')
}

summarize_table_local <- function(input_table,
                                  lm.model = {
                                    trait_value ~ 1 +
                                      group_number +
                                      treatment
                                  },
                                  out.path = 'reports/test.xlsx',
                                  debug = IGNORE_REPORTS) {
  if (!debug) {
    simple_summarize <- input_table %>%
      dplyr::group_by(trait, grouping_gene, group_number, treatment) %>%
      dplyr::summarise(
        n_obs = dplyr::n_distinct(v_t_r),
        mean = mean(trait_value, na.rm = TRUE),
        min = min(trait_value, na.rm = TRUE),
        max = max(trait_value, na.rm = TRUE),
        sd = sd(trait_value, na.rm = TRUE),
        cv_percent = 100 * sd(trait_value, na.rm = TRUE) /
          mean(trait_value, na.rm = TRUE),
      ) %>%
      dplyr::ungroup()

    models_table <- input_table %>% # first linear models
      dplyr::group_by(trait, grouping_gene) %>%
      dplyr::do(model = lm(lm.model, data = .))


    lsmeans_table <- models_table %>%
      dplyr::mutate(lsm_each = list(broom::tidy(emmeans::lsmeans(
        model, c('group_number', 'treatment')
      ))))

    hard_summarize <- lsmeans_table %>% tidyr::unnest(lsm_each)

    signif_table <- input_table %>% # first linear models
      dplyr::group_by(trait, grouping_gene) %>%
      dplyr::do(aov = broom::tidy(aov(lm.model, data = .))) %>%
      tidyr::unnest(aov) %>%
      dplyr::filter(term != 'Residuals') %>%
      dplyr::mutate(p.value = p_stars(p.value)) %>%
      tidyr::pivot_wider(
        names_from = term,
        values_from = p.value,
        names_prefix = 'p_',
        id_cols = c(trait, grouping_gene)
      ) %>%
      dplyr::ungroup()

    hard_summarize %>%
      dplyr::left_join(simple_summarize,
                       by = c('trait', 'grouping_gene',
                              'group_number', 'treatment')) %>%
      dplyr::left_join(signif_table, by = c('trait', 'grouping_gene')) %>%
      dplyr::select(!model) -> final_table

    final_table %>%
      openxlsx::write.xlsx(out.path)

    return(final_table)
  }
}

make_report <- function(table,
                        colnames,
                        title = 'defailt',
                        filename = 'default.html',
                        debug = IGNORE_REPORTS) {
  if (!debug) {
    saveRDS(table, file = '.cache/selected_table.rds')
    saveRDS(colnames, file = '.cache/columns.rds')


    file.copy(from = 'templates/report_timeseries.qmd',
              to = 'report_timeseries.qmd',
              overwrite = TRUE)
    quarto::quarto_render(
      'report_timeseries.qmd',
      output_file = filename,
      execute_params = list('report_title' = title)
    )
    file.remove('report_timeseries.qmd')
  }
}

# import planteye table -------------------------------------------------------
phenospex_file <-
  './data/2022-03-24-Wheat_NO3_#1(b3-6)_20220426_data.zip'

planteye_table <- readr::read_csv(phenospex_file) %>%
  janitor::clean_names(.)

# import unit data --------------------------------------------------------

unit_file_1 <- './data/wheat_NO3_new_handmade.csv'
unit_file_2 <- './data/wheat_NO3_new_translation.csv'

unit_table_1 <- readr::read_csv(unit_file_1)
unit_table_2 <- readr::read_csv(unit_file_2)

unit_table <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = 'V.T.R') %>%
  janitor::clean_names() %>%
  dplyr::rename(repetition = `repeat`)

# import groups table -----------------------------------------------------
group_file <- './data/groups.xlsx'

excel <- readxl::read_xlsx(group_file)

sheetname_vector <- readxl::excel_sheets(group_file) %>%
  unlist()

sheet_list <- sheetname_vector %>%
  purrr::map(\(x) readxl::read_xlsx(group_file, sheet = x))

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

# table merge -------------------------------------------------------------

sym_diff <- function(a, b) {
  sort(setdiff(union(a, b), intersect(a, b)))
}

dim(planteye_table) # 58380    49

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
sym_diff(unit_keys, planteye_keys) # legit empty pots

merged_table <- planteye_table %>%
  dplyr::left_join(unit_table, by = c('unit' = 't_x_y'))

dim(merged_table) # 58380    55

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
sym_diff(unit_keys, group_keys) # sym_difference keys seem legit

merged_table <- merged_table %>%
  dplyr::left_join(groups_table, by = 'var')

dim(merged_table) # 58380    58

merged_table <- merged_table %>%
  dplyr::select(-treatment.x) %>%
  dplyr::rename(treatment = treatment.y) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = 'UTC')) %>%
  dplyr::mutate(
    treatment = as.character(treatment),
    repetition = as.character(repetition),
    plant = as.character(plant)
  )

dim(merged_table) # 58380    57

# aggregation -------------------------------------------------------------
merged_table <- merged_table %>%
  dplyr::mutate(hours_from_start = as.numeric(difftime(timestamp, min(timestamp),
                                           units = 'hours'))
                )

dbscan_cluster <- merged_table %>%
  dplyr::select(hours_from_start) %>%
  dplyr::pull() %>%
  matrix(ncol = 1) %>%
  dbscan::dbscan(., eps = hours_eps) %>%
  .$cluster %>%
  forcats::as_factor(.) %>%
  tibble::as_tibble_col(column_name = 'dbscan_cluster')

merged_table <- merged_table %>%
  dplyr::bind_cols(dbscan_cluster) %>%
  dplyr::select(-hours_from_start)

remove(dbscan_cluster)

(merged_table %>%
  dplyr::mutate(even = as.character(as.integer(dbscan_cluster) %% 2)) %>%
  ggplot2::ggplot(.,
  ggplot2::aes(x = timestamp,
               y = digital_biomass_mm3,
               colour = even)
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_x_datetime()) %>%
  ggplot2::ggsave('reports/clusters_before_cluster_filtering.png', .)

# reshape table -----------------------------------------------------------

string_colnames <- merged_table %>%
  dplyr::select(where(\(x) !is.numeric(x) &
                        !lubridate::is.timepoint(x))) %>%
  colnames()

numeric_colnames <- merged_table %>%
  dplyr::select(where(\(x) is.numeric(x))) %>%
  drop() %>%
  colnames()

merged_table <- merged_table %>%
  dplyr::relocate(all_of(string_colnames), .before = timestamp) %>%
  dplyr::mutate(total_numeric = rowSums(pick(where(is.numeric)),
                                        na.rm = TRUE)) %>%
  dplyr::filter(total_numeric != 0) %>%
  dplyr::select(-total_numeric) # remove rows full of zeroes

merged_table <- merged_table %>%
  tidyr::drop_na(p_ngr5, ppd, rht_b1)

merged_table <- merged_table %>%
  tidyr::pivot_longer(all_of(numeric_colnames),
                      names_to = 'trait',
                      values_to = 'trait_value')

merged_table <- merged_table %>%
  tidyr::pivot_longer(all_of(c('p_ngr5', 'ppd', 'rht_b1')),
                      names_to = 'grouping_gene',
                      values_to = 'group_number')

# remove outlier groups ---------------------------------------------------

outliers_timepoints <- merged_table %>% #delete over-3-sigmas
  dplyr::select(-c(grouping_gene, group_number)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(dbscan_cluster, trait) %>%
  dplyr::summarise(trait_value = (mean(trait_value) - min(trait_value)) /
                                      (max(trait_value) - min(trait_value))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(trait_value = {abs(trait_value - mean(trait_value)) /
        sd(trait_value) > 3}) %>%
  tidyr::replace_na(list(trait_value = TRUE)) %>%
  dplyr::filter(trait_value)

merged_table <- merged_table %>%
  dplyr::anti_join(outliers_timepoints, by = c('dbscan_cluster', 'trait'))

remove(outliers_timepoints)

merged_table <- merged_table %>%  #delete too high stds
  dplyr::group_by(dbscan_cluster, trait) %>%
  dplyr::mutate(std = sd(trait_value)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(trait) %>%
  dplyr::filter((std - mean(std)) / mean(std) <= 3) %>%
  dplyr::select(-std)


(merged_table %>%
    tidyr::pivot_wider(names_from = 'trait',
                       values_from = 'trait_value') %>%
    dplyr::select(-c(grouping_gene, group_number)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(even = as.character(as.integer(dbscan_cluster) %% 2)) %>%
    ggplot2::ggplot(.,
                    ggplot2::aes(x = timestamp,
                                 y = digital_biomass_mm3,
                                 colour = even)
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_x_datetime()) %>%
  ggplot2::ggsave('reports/clusters_after_cluster_filtering.png', .)

# remove old variables -------------------------------------------------------

remove(
  list = c(
    'excel',
    'sheet_list',
    'planteye_table',
    'unit_table_1',
    'unit_table_2',
    'unit_table',
    'groups_table'
  )
)

remove(
  list = c(
    'planteye_keys',
    'unit_keys',
    'group_keys',
    'phenospex_file',
    'unit_file_1',
    'unit_file_2',
    'group_file',
    'sheetname_vector'
  )
)

# in each cluster replace tech repeats by mean ----------------------------
merged_table <- merged_table %>%
  dplyr::group_by(dbscan_cluster) %>%
  dplyr::mutate(timestamp = median(timestamp)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(v_t_r, trait, dbscan_cluster) %>%
  dplyr::mutate(trait_value = median(trait_value)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
# selected time interval, filtration --------------------------------------

selected_table <- merged_table %>%
  dplyr::filter(
    timestamp >= as.POSIXct('2022-04-01 00:00:00', tz = 'UTC')
  ) %>%
  dplyr::filter(
    timestamp <= as.POSIXct('2022-04-01 15:00:00', tz = 'UTC')
  )

remove(merged_table)
# report data from the time interval, raw ---------------------------------

make_report(selected_table,
            numeric_colnames,
            title = 'before normalization',
            filename = 'report_timerseries_before_normalization.html')

summarize_table_local(selected_table,
                      out.path = 'reports/before_normalization.xlsx')

# percentage to logit -----------------------------------------------------

fix_perc_imprecision <- \(x) dplyr::case_when((x >= 0) &
                                                (x <= 1.00) ~ x,
                                              (x < 0) &
                                                (x >= -0.01) ~ 0,
                                              (x > 1) &
                                                (x <= 1.01) ~ 1,
                                              .default = NA)

logit_table <- selected_table %>%
  dplyr::mutate(trait_value = ifelse(
    stringi::stri_detect_fixed(trait,
                               '_percent'),
    stats::qlogis(fix_perc_imprecision(trait_value)),
    trait_value
  )) %>%
  dplyr::mutate(trait = stringi::stri_replace_all_fixed(trait,
                                                        pattern = 'percent',
                                                        replacement = 'logit')) # logit transform for all percent data, fixing percents from range -0.01 to 1.01

logit_table %>%
  dplyr::filter(stringi::stri_detect_fixed(trait, 'logit')) %>%
  dplyr::summarise(min = min(trait_value),
                   max = max(trait_value)) # logits contain -Inf

# will replace -Inf with closes negative value

minus_inf_replacement <- logit_table %>%
  dplyr::filter(stringi::stri_detect_fixed(trait, 'logit')) %>%
  dplyr::filter(trait_value > -Inf) %>%
  dplyr::summarise(min = min(trait_value)) %>%
  dplyr::pull() %>%
  floor()

logit_table <- logit_table %>%
  dplyr::mutate(trait_value = ifelse(trait_value == -Inf,
                                     minus_inf_replacement,
                                     trait_value))

logit_numeric_colnames <- logit_table %>%
  dplyr::select(trait) %>%
  dplyr::distinct(trait) %>%
  dplyr::pull()

make_report(logit_table,
            logit_numeric_colnames,
            title = 'logit',
            filename = 'report_timerseries_logit.html')

summarize_table_local(logit_table,
                      out.path = 'reports/logit.xlsx')

# distribution in each subgroup -------------------------------------------

logit_gaus_table <-
  logit_table %>% # normalize for each trait
  dplyr::group_by(trait) %>%
  dplyr::mutate(trait_value = LambertW::Gaussianize(trait_value,
                                                    type = 'hh',
                                                    method = 'MLE')) %>%
  dplyr::ungroup()

# MODELS -----------------------------------------------------------------------

logit_gaus_table <- logit_gaus_table %>%
  dplyr::group_by(trait, grouping_gene) %>%
  dplyr::mutate(degenerate = (dplyr::n_distinct(group_number) - 1) *
                  (dplyr::n_distinct(treatment) - 1) == 0) %>%
  dplyr::filter(degenerate == FALSE) %>%
  dplyr::select(-degenerate) %>%
  dplyr::ungroup()

make_report(logit_gaus_table,
            logit_numeric_colnames,
            title = 'after normalization',
            filename = 'report_timerseries_after_normalization.html')

summarize_table_local(logit_gaus_table,
                      out.path = 'reports/after_normalization.xlsx')

# remove redundant tables -------------------------------------------------
remove(list = c(
  'selected_table',
  'logit_table'
))
remove(list = c('minus_inf_replacement', 'numeric_colnames'))
# violins -----------------------------------------------------------------

nested_table <- logit_gaus_table %>%
  dplyr::nest_by(trait, grouping_gene, .key = 'subdf')

ANOVA_table <- nested_table %>%
  dplyr::do(ANOVA = aov(trait_value ~
                          group_number * treatment +
                          1 / var, data = .$subdf))

THSD_table <- ANOVA_table %>%
  dplyr::do(THSD = TukeyHSD(.$ANOVA))

THSD_table <- THSD_table %>%
  as.list %>%
  purrr::modify_depth(5, \(x)
                      tidyr::replace_na(x, replace = 0),
                      .ragged = TRUE) %>%
  tibble::as_tibble() #fill na-ridden rows with 0 diff 0 p-val

CL_table <- ANOVA_table %>%
  dplyr::bind_cols(THSD_table) %>%
  dplyr::do(CL = multcompView::multcompLetters4(.$ANOVA, .$THSD))

letter_table <- nested_table %>%
  dplyr::select(-subdf) %>%
  dplyr::bind_cols(CL_table) %>%
  dplyr::do(tukey_groupings = tidy_CL_cell(.$CL))

tukeys_groups_table <- nested_table %>%
  dplyr::select(-subdf) %>%
  dplyr::bind_cols(letter_table) %>%
  tidyr::unnest(tukey_groupings)

remove(list = c(
  'nested_table',
  'ANOVA_table',
  'THSD_table',
  'CL_table',
  'letter_table'
))

printable_table <- logit_gaus_table %>%
  dplyr::mutate(group_number_treatment = paste(group_number,
                                               treatment,
                                               sep = ':')) %>%
  tidyr::pivot_longer(
    c(group_number, treatment, group_number_treatment),
    names_to = 'tukey_grouping',
    values_to = 'tukey_group'
  ) %>%
  dplyr::left_join(
    tukeys_groups_table,
    by = c('trait', 'grouping_gene',
           'tukey_grouping', 'tukey_group')
  )

gc()
if (!IGNORE_REPORTS) {
  named_plots <- printable_table %>%
    dplyr::group_by(trait, grouping_gene, tukey_grouping) %>%
    dplyr::arrange(tukey_group) %>%
    dplyr::group_map(
      ~ list(
        paste(.y[[1, 'trait']],
              .y[[1, 'grouping_gene']],
              .y[[1, 'tukey_grouping']],
              sep = ', '),
        ggplot2::ggplot(
          .,
          ggplot2::aes(y = trait_value,
                       x = tukey_group,
                       color = tukey_letter)
        ) +
          ggplot2::geom_violin() +
          ggplot2::geom_boxplot(width = 0.1) +
          ggplot2::ggtitle(paste(.y[[1, 'trait']],
                                 .y[[1, 'grouping_gene']],
                                 .y[[1, 'tukey_grouping']],
                                 sep = ', ')) +
          ggplot2::theme_minimal()
      )
    )

  saveRDS(named_plots, file = '.cache/named_plots.rds')

  file.copy(from = 'templates/report_violins.qmd',
            to = 'report_violins.qmd',
            overwrite = TRUE)
  quarto::quarto_render(
    'report_violins.qmd',
    output_file = 'violins_gaus.html',
    execute_params = list('report_title' = 'Violins')
  )
  file.remove('report_violins.qmd')
}
