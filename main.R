
# install libraries -------------------------------------------------------
install.packages('tidyverse')
install.packages('janitor')
install.packages('readxl')
install.packages('stringi')
install.packages('emmeans')
install.packages('lme4')


# load libraries ----------------------------------------------------------
library('tidyverse')
library('janitor')
library('readxl')
library('stringi')
library('emmeans')
library('lme4')

# import planteye table -------------------------------------------------------
table_file <- './data/2022-03-24-Wheat_NO3_#1(b3-6)_20220426_data.zip'

planteye_table <- readr::read_csv(table_file) %>% 
  janitor::clean_names(.)

# import unit data --------------------------------------------------------

unit_file_1 <- './data/wheat_NO3_new_handmade.csv'
unit_file_2 <- './data/wheat_NO3_new_translation.csv'

unit_table_1 <- readr::read_csv(unit_file_1)
unit_table_2 <- readr::read_csv(unit_file_2)

unit_table <- unit_table_1 %>% 
  dplyr::left_join(unit_table_2, by='V.T.R') %>% 
  janitor::clean_names() %>% 
  dplyr::rename(repetition=`repeat`)

# import groups table -----------------------------------------------------
group_file <- './data/groups.xlsx'

excel <- readxl::read_xlsx(group_file)

sheetname_vector <- readxl::excel_sheets(group_file) %>% 
  unlist

sheet_list <- sheetname_vector %>% 
  purrr::map(\(x) readxl::read_xlsx(group_file, sheet=x)) 

groups_table <-
  purrr::map2_dfr(sheetname_vector,
                  sheet_list,
                  \(naming, subtable) subtable %>% 
                    dplyr::mutate(source = naming))

groups_table <- groups_table %>% 
  pivot_longer(contains('Вариант'), 
               names_to = 'group_number', 
               values_to = 'var') %>% 
  dplyr::select(!contains('...')) %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(
    group_number = stringi::stri_replace_all_regex(group_number, '[^\\d]','')
    ) %>% 
  dplyr::rename(group_principle=source) %>% 
  pivot_wider(names_from = group_principle,
              values_from = group_number) %>% 
  janitor::clean_names(.)


# table merge -------------------------------------------------------------

sym_diff <- function(a,b) sort(setdiff(union(a,b), intersect(a,b)))

dim(planteye_table) #58380    49

planteye_keys <- planteye_table %>% 
  select(unit) %>% 
  distinct() %>% 
  pull %>% 
  sort
unit_keys <- unit_table %>% 
  select(t_x_y) %>% 
  distinct() %>% 
  pull %>% 
  sort
sym_diff(unit_keys, planteye_keys) #legit empty pots

merged_table <- planteye_table %>% 
  dplyr::left_join(unit_table, by=c('unit'='t_x_y'))

dim(merged_table) #58380    55

unit_keys <- unit_table %>% select(var) %>% distinct() %>% pull %>% sort
group_keys <- groups_table %>% select(var) %>% distinct() %>% pull %>% sort
sym_diff(unit_keys, group_keys) #sym_difference keys seem legit


merged_table <- merged_table %>% #58380    58
  dplyr::left_join(groups_table, by='var')

dim(merged_table) #58380    58

merged_table %>% 
  select(-'treatment.x') %>% 
  rename(treatment = treatment.y)

# selected time interval, filtration --------------------------------------

selected_table <- merged_table %>% 
  filter(timestamp >= as.Date("2022-04-01"), 
         timestamp < as.Date("2022-04-02"))

selected_table <- selected_table %>% 
  select(-treatment.x) %>% 
  rename(treatment=treatment.y) %>% 
  mutate(treatment=as.character(treatment),
         repetition=as.character(repetition),
         plant=as.character(plant))

string_colnames <- selected_table %>% 
  select(where(\(x) !is.numeric(x) & !is.timepoint(x))) %>% 
  colnames

numeric_colnames <- selected_table %>% 
  select(where(\(x) is.numeric(x))) %>%
  drop() %>% 
  colnames

selected_table <- selected_table %>% 
  relocate(all_of(string_colnames), .before=timestamp)

selected_table <- selected_table %>% 
  mutate(total_numeric = rowSums(pick(where(is.numeric)), na.rm=TRUE)) %>% 
  filter(total_numeric != 0)

selected_table <- selected_table %>% 
  drop_na(p_ngr5, ppd, rht_b1) #ONLY GROUPED VARS, 912 × 58

selected_table <- selected_table %>% 
  pivot_longer(all_of(numeric_colnames), 
               names_to = 'trait', 
               values_to= 'trait_value')

selected_table <- selected_table %>% 
  pivot_longer(all_of(c('p_ngr5', 'ppd', 'rht_b1')), 
               names_to = 'grouping_gene', 
               values_to= 'group_number')

simple_summarize <- selected_table %>% 
  group_by(trait, grouping_gene, group_number, treatment) %>% 
  summarise(n_obs = n_distinct(v_t_r),
            mean = mean(trait_value, na.rm=TRUE),
            min = min(trait_value, na.rm=TRUE),
            max = max(trait_value, na.rm=TRUE),
            sd = sd(trait_value, na.rm=TRUE),
            cv_percent = 100*mean(trait_value, na.rm=TRUE)/sd(trait_value, na.rm=TRUE),
            )

models <- selected_table %>% #first linear models
  group_by(trait, grouping_gene) %>%
  do(model = lm(trait_value ~ 0 + group_number + treatment, data=.))