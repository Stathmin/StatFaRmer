# init --------------------------------------------------------------------
renv::activate()
renv::hydrate(prompt = FALSE)
here::i_am("README.md")
`%>%` <- magrittr::`%>%`

IGNORE_REPORTS = FALSE

# functions ---------------------------------------------------------------
summarize_table_local <- function(input_table,
                                  lm.model = {
                                    trait_value ~ 0 +
                                      group_number +
                                      treatment
                                  },
                                  out.path = "reports/test.xlsx",
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
      )
    
    models_table <- input_table %>% # first linear models
      dplyr::group_by(trait, grouping_gene) %>%
      dplyr::do(model = lm(lm.model, data = .))
    
    
    lsmeans_table <- models_table %>%
      dplyr::mutate(lsm_each = list(broom::tidy(emmeans::lsmeans(
        model, c("group_number", "treatment")
      ))))
    
    hard_summarize <- lsmeans_table %>% tidyr::unnest(lsm_each)
    
    signif_table <- input_table %>% # first linear models
      dplyr::group_by(trait, grouping_gene) %>%
      dplyr::do(aov = broom::tidy(aov(lm.model, data = .))) %>%
      tidyr::unnest(aov) %>%
      dplyr::filter(term != "Residuals") %>%
      dplyr::mutate(
        p.value = dplyr::case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          p.value < 0.1 ~ ".",
          .default = " "
        )
      ) %>%
      tidyr::pivot_wider(
        names_from = term,
        values_from = p.value,
        names_prefix = "p_",
        id_cols = c(trait, grouping_gene)
      )
    
    hard_summarize %>%
      dplyr::left_join(simple_summarize,
                       by = c("trait", "grouping_gene",
                              "group_number", "treatment")) %>%
      dplyr::left_join(signif_table, by = c("trait", "grouping_gene")) %>%
      dplyr::select(!model) -> final_table
    
    final_table %>%
      openxlsx::write.xlsx(out.path)
    
    return(final_table)
  }
}

make_report <- function(table,
                        colnames,
                        title = "defailt",
                        filename = "default.html",
                        debug = IGNORE_REPORTS) {
  if (!debug) {
    saveRDS(table, file = ".cache/selected_table.rds")
    saveRDS(colnames, file = ".cache/columns.rds")
    
    
    file.copy(from = "templates/report_timeseries.qmd",
              to = "report_timeseries.qmd",
              overwrite = TRUE)
    quarto::quarto_render(
      "report_timeseries.qmd",
      output_file = filename,
      execute_params = list("report_title" = title)
    )
    file.remove("report_timeseries.qmd")
  }
}

# import planteye table -------------------------------------------------------
table_file <-
  "./data/2022-03-24-Wheat_NO3_#1(b3-6)_20220426_data.zip"

planteye_table <- readr::read_csv(table_file) %>%
  janitor::clean_names(.)

# import unit data --------------------------------------------------------

unit_file_1 <- "./data/wheat_NO3_new_handmade.csv"
unit_file_2 <- "./data/wheat_NO3_new_translation.csv"

unit_table_1 <- readr::read_csv(unit_file_1)
unit_table_2 <- readr::read_csv(unit_file_2)

unit_table <- unit_table_1 %>%
  dplyr::left_join(unit_table_2, by = "V.T.R") %>%
  janitor::clean_names() %>%
  dplyr::rename(repetition = `repeat`)

# import groups table -----------------------------------------------------
group_file <- "./data/groups.xlsx"

excel <- readxl::read_xlsx(group_file)

sheetname_vector <- readxl::excel_sheets(group_file) %>%
  unlist()

sheet_list <- sheetname_vector %>%
  purrr::map(\(x) readxl::read_xlsx(group_file, sheet = x))

groups_table <-
  purrr::map2_dfr(
    sheetname_vector,
    sheet_list,
    \(naming, subtable) subtable %>%
      dplyr::mutate(source = naming)
  )

groups_table <- groups_table %>%
  tidyr::pivot_longer(contains("Вариант"),
    names_to = "group_number",
    values_to = "var"
  ) %>%
  dplyr::select(!contains("...")) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(group_number = stringi::stri_replace_all_regex(group_number, "[^\\d]", "")) %>%
  dplyr::rename(group_principle = source) %>%
  tidyr::pivot_wider(
    names_from = group_principle,
    values_from = group_number
  ) %>%
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
  dplyr::left_join(unit_table, by = c("unit" = "t_x_y"))

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
  dplyr::left_join(groups_table, by = "var")

dim(merged_table) # 58380    58

merged_table <- merged_table %>%
  dplyr::select(-treatment.x) %>%
  dplyr::rename(treatment = treatment.y) %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) %>%
  dplyr::mutate(
    treatment = as.character(treatment),
    repetition = as.character(repetition),
    plant = as.character(plant)
  )

dim(merged_table) # 58380    57


# selected time interval, filtration --------------------------------------

selected_table <- merged_table %>%
  dplyr::filter(
    timestamp >= as.POSIXct("2022-04-01 00:00:00", tz = "UTC"),
    timestamp <= as.POSIXct("2022-04-01 15:00:00", tz = "UTC")
  )

string_colnames <- selected_table %>%
  dplyr::select(where(\(x) !is.numeric(x) &
    !lubridate::is.timepoint(x))) %>%
  colnames()

numeric_colnames <- selected_table %>%
  dplyr::select(where(\(x) is.numeric(x))) %>%
  drop() %>%
  colnames()

selected_table <- selected_table %>%
  dplyr::relocate(all_of(string_colnames), .before = timestamp) %>%
  dplyr::mutate(total_numeric = rowSums(pick(where(is.numeric)),
    na.rm = TRUE
  )) %>%
  dplyr::filter(total_numeric != 0) %>%
  dplyr::select(-total_numeric) # remove rows full of zeroes

selected_table <- selected_table %>%
  tidyr::drop_na(p_ngr5, ppd, rht_b1) # KEEP ONLY PLANTS WITH KNOWN GT, 456  57


# pivot selected table ---------------------------------------------------------

selected_table <- selected_table %>%
  tidyr::pivot_longer(all_of(numeric_colnames),
    names_to = "trait",
    values_to = "trait_value"
  )

selected_table <- selected_table %>%
  tidyr::pivot_longer(all_of(c("p_ngr5", "ppd", "rht_b1")),
    names_to = "grouping_gene",
    values_to = "group_number"
  )

# report data from the time interval, raw ---------------------------------

make_report(selected_table,
  numeric_colnames,
  title = "before normalization",
  filename = "report_timerseries_before_normalization.html"
)

summarize_table_local(selected_table,
  out.path = "reports/before_normalization.xlsx"
)


# percentage to logit -----------------------------------------------------

fix_perc_imprecision <- \(x) dplyr::case_when(
  (x >= 0) &
    (x <= 1.00) ~ x,
  (x < 0) & (x >= -0.01) ~ 0,
  (x > 1) & (x <= 1.01) ~ 1,
  .default = NA
)

logit_table <- selected_table %>%
  dplyr::mutate(trait_value = ifelse(
    stringi::stri_detect_fixed(
      trait,
      "_percent"
    ),
    stats::qlogis(fix_perc_imprecision(trait_value)),
    trait_value
  )) %>%
  dplyr::mutate(trait = stringi::stri_replace_all_fixed(trait,
    pattern = "percent",
    replacement = "logit"
  )) # logit transform for all percent data, fixing percents from range -0.01 to 1.01

logit_table %>%
  dplyr::filter(stringi::stri_detect_fixed(trait, "logit")) %>%
  dplyr::summarise(
    min = min(trait_value),
    max = max(trait_value)
  ) # logits contain -Inf

# will replace -Inf with closes negative value

minus_inf_replacement <- logit_table %>%
  dplyr::filter(stringi::stri_detect_fixed(trait, "logit")) %>%
  dplyr::filter(trait_value > -Inf) %>%
  dplyr::summarise(min = min(trait_value)) %>%
  dplyr::pull() %>%
  floor()

logit_table <- logit_table %>%
  dplyr::mutate(trait_value = ifelse(trait_value == -Inf,
    minus_inf_replacement,
    trait_value
  ))

logit_numeric_colnames <- logit_table %>%
  dplyr::select(trait) %>%
  dplyr::distinct() %>%
  dplyr::pull()


make_report(logit_table,
  logit_numeric_colnames,
  title = "logit",
  filename = "report_timerseries_logit.html"
)

summarize_table_local(logit_table,
  out.path = "reports/logit.xlsx"
)

# distribution in each subgroup -------------------------------------------

logit_3sigma_table <- logit_table %>%
  dplyr::group_by(var, treatment, trait) %>%
  dplyr::mutate(z_score = abs(trait_value - mean(trait_value)) / sd(trait_value), ) %>%
  dplyr::filter(z_score <= 3)

logit_gaus_table <- logit_3sigma_table %>% # normalize for each trait
  dplyr::group_by(trait) %>%
  dplyr::mutate(trait_value = LambertW::Gaussianize(trait_value))

# MODELS -----------------------------------------------------------------------

logit_gaus_table <- logit_gaus_table %>%
  dplyr::group_by(trait, grouping_gene) %>%
  dplyr::mutate(degenerate = (dplyr::n_distinct(group_number) - 1) *
    (dplyr::n_distinct(treatment) - 1) == 0) %>%
  dplyr::filter(degenerate == FALSE) %>%
  dplyr::select(-degenerate)

make_report(logit_gaus_table,
  logit_numeric_colnames,
  title = "after normalization",
  filename = "report_timerseries_after_normalization.html"
)

summarize_table_local(logit_gaus_table,
  out.path = "reports/after_normalization.xlsx"
)


# violins -----------------------------------------------------------------
named_plots <- logit_gaus_table %>% 
  dplyr::group_by(trait, grouping_gene) %>% 
  dplyr::arrange(group_number, treatment) %>% 
  dplyr::mutate(gr_tr = paste(group_number, treatment, sep = '_')) %>% 
  dplyr::group_map(
    ~ list(paste(.y[[1,'trait']],
                 .y[[1,'grouping_gene']],
                 sep = ', '),
      ggplot2::ggplot(., ggplot2::aes(y = trait_value, 
                                 x = gr_tr, 
                                 color = group_number)) + 
      ggplot2::geom_violin() +
      ggplot2::geom_boxplot(width = 0.1) +
      ggplot2::ggtitle(paste(.y[[1,'trait']],
                             .y[[1,'grouping_gene']],
                             sep = ', ')) +
      ggplot2::theme_minimal())
  )
  
saveRDS(named_plots, file = ".cache/named_plots.rds")


file.copy(from = "templates/report_violins.qmd",
          to = "report_violins.qmd",
          overwrite = TRUE)
quarto::quarto_render(
  "report_violins.qmd",
  output_file = 'violins_gaus.html',
  execute_params = list("report_title" = 'Violins')
)
file.remove("report_violins.qmd")
