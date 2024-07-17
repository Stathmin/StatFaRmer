library('shiny')
library('tidyverse')
library('stringr')
library('ggplot2')
library('multcompView')
library('broom')
library('flextable')
library('moments')
library('bslib')
library('viridis')
library('svglite')

options(shiny.reactlog = TRUE)
set.seed(42)

formulate <- function(RHS, LHF) {
  if (length(LHF) == 0) {
    LHS = '1'
  } else if ((length(LHF) == 1)) {
    LHS = str_interp('1 + ${LHF}')
  } else {
    LHS = str_interp('1 + (${paste(LHF, collapse=" + ")})^2')
  }
  return(paste(RHS, LHS, sep = " ~ "))
}

generate_label_df <- function(TUKEY, variable) {
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][, 4]
  Tukey.labels <-
    data.frame(multcompLetters(Tukey.levels)['Letters'])

  #Keep labels in the same order as in the boxplot :
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment) , ]
  colnames(Tukey.labels) = c('letter', 'group')
  return(Tukey.labels)
}


merged_table <-
  readRDS(stringr::str_interp('merged_table.rds'))
vector_of_groups <-
  readRDS(stringr::str_interp('vector_of_groups.rds'))
merged_table[sapply(merged_table, is.infinite)] <- NA
merged_table <- merged_table %>%
  arrange(timestamp) %>%
  group_by(dbscan_cluster) %>%
  mutate(timestamp_group = mean(timestamp)) %>%
  ungroup()

get_unique <- function(table, string) {
  table %>%
    select(!!string) %>%
    distinct() %>%
    pull %>%
    sort
}

grouping_factors <- merged_table %>%
  select(!is.numeric, -timestamp) %>%
  colnames %>% sort

treatments <- get_unique(merged_table, 'treatment')
cultivars <- get_unique(merged_table, 'cultivar')
timestamp_groups <- get_unique(merged_table, 'timestamp_group')
named_timestamp_groups <- timestamp_groups
names(named_timestamp_groups) <-
  timestamp_groups %>% format(format = '%m-%d %H:%M')

out_variables <-
  merged_table %>% select(where(is.numeric)) %>% colnames() %>% sort()



# Define UI for application that draws a histogram
ui <- fluidPage(# Application title
  titlePanel("StatFaRmer"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput(
        'grouping_factors',
        'Grouping factors:',
        grouping_factors[grouping_factors != "timestamp_group"],
        multiple = TRUE,
        selected = c('treatment', 'dbscan_clusters')
      ),
      selectInput(
        'tukey_group',
        'Tukey group:',
        grouping_factors[grouping_factors != "timestamp_group"],
        multiple = TRUE,
        selected = c('treatment', 'dbscan_clusters')
      ),
      selectInput(
        'gene_grouping',
        'Grouping gene:',
        choices = vector_of_groups,
        selected = {
          vector_of_groups %>% tail(1)
        }
      ),
      selectizeInput(
        'gene_groups',
        'Gene groups:',
        choices = get_unique(merged_table, {
          vector_of_groups %>% tail(1)
        }),
        multiple = TRUE,
        selected = get_unique(merged_table, {
          vector_of_groups %>% tail(1)
        })
      ),
      selectInput(
        'treatments',
        'Selected treatments:',
        choices = treatments,
        multiple = TRUE,
        selected = treatments
      ),
      selectInput(
        'cultivars',
        'Selected cultivars:',
        choices = cultivars,
        multiple = TRUE,
        selected = cultivars
      ),
      shinyWidgets::pickerInput(
        inputId = "timestamp_groups",
        label = "Selected time clusters:",
        choices = named_timestamp_groups,
        selected = named_timestamp_groups[c(1, round(length(named_timestamp_groups) /
                                                       2), length(named_timestamp_groups))] %>%
          as.character(),
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, size = 10),
        multiple = TRUE
      ),
      selectInput(
        'out_variables',
        'Selected trait:',
        choices = out_variables,
        multiple = FALSE,
        selected = out_variables[1]
      ),
      textInput(
        'facet_formula',
        'Facet formula:',
        value = 'treatment ~ timestamp_group',
        placeholder = 'treatment ~ timestamp_group'
      ),
      checkboxInput('timeseries_plot', ': plot timeseries with GAM', FALSE),

      actionButton("submit", "Submit"),
      downloadButton("savePlot", "Save Plot as SVG")
    ),

    mainPanel(
      uiOutput("formula"),
      plotOutput("distPlot", height = "800px", width = '1200px'),
      tabsetPanel(
        tabPanel(
          "Raw Table",
          uiOutput("RAW"),
          verbatimTextOutput("raw"),
          DT::DTOutput("raw_flex"),
          downloadButton("raw_flex_downloadData", "Download Full Results")
        ),
        tabPanel(
          "Descriptive",
          uiOutput("DESCvarNames"),
          verbatimTextOutput("DescriptiveTable"),
          DT::DTOutput("descriptiveTable_flex"),
          downloadButton("descriptiveTable_flex_downloadData", "Download Full Results")
        ),
        tabPanel(
          "ANOVA",
          uiOutput("ANOVAvarNames"),
          verbatimTextOutput("anovaTable"),
          DT::DTOutput("anovaTable_flex"),
          downloadButton("anovaTable_flex_downloadData", "Download Full Results")
        ),
        tabPanel(
          "Tukey",
          uiOutput("TUKEYvarNames"),
          verbatimTextOutput("tukeyTable"),
          DT::DTOutput("tukeyTable_flex"),
          downloadButton("tukeyTable_flex_downloadData", "Download Full Results")
        ),
        tabPanel(
          "Group Letters",
          uiOutput("TUKEYLvarNames"),
          verbatimTextOutput("tukeyLetters"),
          DT::DTOutput("tukeyLetters_flex"),
          downloadButton("tukeyLetters_flex_downloadData", "Download Full Results")
        )
      )
    )
  ))

# Define server logic required to draw a histogram


server <- function(input, output, session) {
  set.seed(42)

  combined_inputs <- reactiveValues()

  create_initial_combined_inputs <- reactive({
    list(
      gene_grouping = isolate(input$gene_grouping),
      gene_groups = isolate(input$gene_groups),
      treatments = isolate(input$treatments),
      cultivars = isolate(input$cultivars),
      timestamp_groups = isolate(unname(input$timestamp_groups)),
      out_variables = isolate(input$out_variables),
      facet_formula = isolate(as.formula(input$facet_formula)),
      grouping_factors = isolate(sort(input$grouping_factors)),
      tukey_group = isolate(sort(input$tukey_group)),
      timeseries_plot = isolate(input$timeseries_plot)
    )
  })

  observe({
    initial_combined_inputs <- create_initial_combined_inputs()

    combined_inputs$gene_grouping <- initial_combined_inputs$gene_grouping
    combined_inputs$gene_groups <- initial_combined_inputs$gene_groups
    combined_inputs$treatments <- initial_combined_inputs$treatments
    combined_inputs$cultivars <- initial_combined_inputs$cultivars
    combined_inputs$timestamp_groups <- initial_combined_inputs$timestamp_groups
    combined_inputs$out_variables <- initial_combined_inputs$out_variables
    combined_inputs$facet_formula <- initial_combined_inputs$facet_formula
    combined_inputs$grouping_factors <- initial_combined_inputs$grouping_factors
    combined_inputs$tukey_group <- initial_combined_inputs$tukey_group
    combined_inputs$timeseries_plot <- initial_combined_inputs$timeseries_plot
  }, priority = 50)

  observeEvent(input$submit, {
    combined_inputs$gene_grouping = isolate(input$gene_grouping)
    combined_inputs$gene_groups = isolate(input$gene_groups)
    combined_inputs$treatments = isolate(input$treatments)
    combined_inputs$cultivars = isolate(input$cultivars)
    combined_inputs$timestamp_groups = isolate(unname(input$timestamp_groups))
    combined_inputs$out_variables = isolate(input$out_variables)
    combined_inputs$facet_formula = isolate(as.formula(input$facet_formula))
    combined_inputs$grouping_factors = isolate(sort(input$grouping_factors))
    combined_inputs$tukey_group = isolate(sort(input$tukey_group))
    combined_inputs$timeseries_plot = isolate(input$timeseries_plot)
  }, priority = 50)

  observeEvent(input$gene_grouping,
               {
                 {
                   unique_groups <- get_unique(merged_table, input$gene_grouping)

                   if (!setequal(unique_groups, input$gene_groups)) {
                     isolate(
                       updateSelectizeInput(
                         session,
                         'gene_groups',
                         choices = c(),
                         selected = c(),
                         server = TRUE
                       )
                     )
                     updateSelectizeInput(
                       session,
                       'gene_groups',
                       choices = unique_groups,
                       selected = unique_groups,
                       server = TRUE
                     )
                   }
                 }
               },
               ignoreInit = FALSE,
               priority = 100)

  observeEvent(
    c(input$gene_groups, input$gene_grouping),
    {
      {
        unique_cultivars <- merged_table %>%
          filter(!!sym(input$gene_grouping) %in% input$gene_groups) %>%
          select('cultivar') %>%
          distinct() %>%
          arrange() %>%
          pull()

        if (!setequal(unique_cultivars, input$cultivars)) {
          isolate(
            updateSelectizeInput(
              session,
              'cultivars',
              choices = c(),
              selected = c(),
              server = TRUE
            )
          )
          updateSelectizeInput(
            session,
            'cultivars',
            choices = unique_cultivars,
            selected = unique_cultivars,
            server = TRUE
          )
        }
      }
    },
    ignoreInit = FALSE,
    priority = 75
  )

  observeEvent(
    c(
      combined_inputs$treatments,
      combined_inputs$cultivars,
      combined_inputs$timestamp_groups,
      combined_inputs$out_variables,
      combined_inputs$facet_formula,
      combined_inputs$grouping_factors,
      combined_inputs$tukey_group
    ),
    {
        {
        local_model_d <- {
          formulate(combined_inputs$out_variables,
                    combined_inputs$grouping_factors)
        }
        logger::log_info(paste0('local_model_d:', local_model_d))

        most_interactive <- {
          as.formula(local_model_d) %>%
            terms() %>%
            labels() %>%
            {
              if (length(.) == 0) {
                FALSE
              } else {
                .[str_count(., ':') == max(str_count(., ':'))][1] %>%
                  str_split_1(':') %>%
                  unlist()
              }
            }
        }
        logger::log_info(paste0('most_interactive:', paste0(most_interactive, collapse = ', ')))

        initial_table <- {
          merged_table %>%
            filter(
              treatment %in% combined_inputs$treatments,
              cultivar %in% combined_inputs$cultivars,
              as.character(timestamp_group) %in% as.character(combined_inputs$timestamp_groups),
              !!sym(combined_inputs$gene_grouping) %in% combined_inputs$gene_groups
            ) %>%
            arrange(timestamp_group) %>%
            mutate(timestamp_group = as.factor(timestamp_group)) %>%
            {
              if (!is.null(combined_inputs$tukey_group)) {
                unite(
                  data = .,
                  col = group,
                  !!(sort(combined_inputs$tukey_group)),
                  sep = ":",
                  remove = FALSE
                )
              } else {
                mutate(., group = 'all')
              }
            } %>%
            select(-c(where(is.numeric), -!!combined_inputs$out_variables))
        }
        logger::log_info(paste0('length_initial_table:', nrow(initial_table)))
        }#fast parse

        {
          local_descriptive <- {
            initial_table %>%
              filter(!is.na(!!sym(
                combined_inputs$out_variables
              ))) %>%
              arrange(!!sym(combined_inputs$out_variables)) %>%
              {
                if (!is.null(combined_inputs$tukey_group)) {
                  unite(
                    data = .,
                    col = group,
                    !!(sort(combined_inputs$tukey_group)),
                    sep = ":",
                    remove = FALSE
                  )
                } else {
                  mutate(., group = 'all')
                }
              } %>%
              select(group, !!combined_inputs$out_variables) %>%
              group_by(group) %>%
              summarise(
                n = dplyr::n(),
                median = median(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                mean = mean(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                cv_perc = 100 * sd(!!rlang::sym(
                  combined_inputs$out_variables
                )) / median(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                min = min(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                max = max(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                skewness = moments::skewness(!!rlang::sym(
                  combined_inputs$out_variables
                )),
                kurtosis = moments::kurtosis(!!rlang::sym(
                  combined_inputs$out_variables
                )),
              ) %>%
              ungroup() %>%
              arrange(mean)
          }
          logger::log_info(paste0('length_local_descriptive:', nrow(local_descriptive)))

          local_anova <- {
            aov(as.formula(local_model_d), data = initial_table)
          }
          tukey_needed <- {
            as.formula(local_model_d) %>%
              all.vars() %>%
              length() > 1
          }
          logger::log_info(paste0('tukey_needed:', tukey_needed))

          letters_needed <- {
            tukey_needed &
              length(vecsets::vsetdiff(
                sort(combined_inputs$tukey_group),
                sort(combined_inputs$grouping_factors)
              )) == 0
          }
          logger::log_info(paste0('letters_needed:', letters_needed))

          local_tukey <- if (tukey_needed) {
            TukeyHSD(local_anova, ordered = TRUE) %>%
              purrr::modify_depth(5, ~ replace_na(.x, replace = 0), .ragged = TRUE)
          } else {
            NULL
          }

          local_names <- if (letters_needed) {
            multcompLetters4(local_anova, local_tukey) %>%
              .[[paste(sort(combined_inputs$tukey_group),
                       sep = ':',
                       collapse = ':')]] %>%
              as.data.frame.list() %>%
              as_tibble(rownames = 'group') %>%
              select('group', 'Letters') %>%
              rename(letter = Letters)
          } else {
            NULL
          }

          local_table <- if (letters_needed) {
            initial_table %>%
              filter(!is.na(!!combined_inputs$out_variables)) %>%
              left_join(local_names, by = c('group' = 'group'))
          } else {
            initial_table %>%
              mutate(letter = '-')
          }
          logger::log_info(paste0('local_table:', nrow(local_table)))
        }#anova-tukey-letters-calc

        {
          output$formula <- renderUI({
            HTML(as.character(
              div(style = "text-align: center; font-weight: bold;", local_model_d)
            ))
          })

          render_local_table <- local_table %>%
            mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

          output$raw_flex <- DT::renderDT(
            DT::datatable(
              render_local_table,
              options = list(
                server = TRUE,
                filter = "top",
                selection = "multiple",
                searching = TRUE,
                fixedColumns = TRUE,
                autoWidth = TRUE,
                ordering = TRUE,
                dom = 'tBp',
                buttons = list(
                  'excel', 'csv')
              ),
              class = "display",
              extensions = 'Buttons'
            ),
            server = TRUE
          )


        output$raw_flex_downloadData <- downloadHandler(
          filename = function() {
            paste("Raw_", Sys.Date(), ".xlsx", sep = "")
          },
          content = function(file) {
            # Write the full dataset to an Excel file
            writexl::write_xlsx(render_local_table, file)
          }
        )

        render_descriptive <- local_descriptive %>%
          mutate(across(where(is.numeric), ~ round(., digits = 2)))

        output$descriptiveTable_flex <- DT::renderDT(
          DT::datatable(
            render_descriptive,
            options = list(
              server = TRUE,
              filter = "top",
              selection = "multiple",
              searching = TRUE,
              fixedColumns = TRUE,
              autoWidth = TRUE,
              ordering = TRUE,
              dom = 'tBp',
              buttons = list(
                'excel', 'csv')
            ),
            class = "display",
            extensions = 'Buttons'
          ),
          server = TRUE
        )

        output$descriptiveTable_flex_downloadData <- downloadHandler(
          filename = function() {
            paste("Descriptive_", Sys.Date(), ".xlsx", sep = "")
          },
          content = function(file) {
            # Write the full dataset to an Excel file
            writexl::write_xlsx(render_descriptive, file)
          }
        )

        render_local_anova <- local_anova %>% tidy() %>%
          mutate(
            sig = case_when(
              p.value <= 0.001 ~ "***",
              p.value <= 0.01 ~ "**",
              p.value <= 0.05 ~ "*",
              p.value <= 0.1 ~ ".",
              TRUE ~ ""
            )
          ) %>%
          mutate(across(where(is.numeric), ~ round(., digits = 2)))

          output$anovaTable_flex <- {DT::renderDT(
            DT::datatable(
              render_local_anova,
              options = list(
                server = TRUE,
                filter = "top",
                selection = "multiple",
                searching = TRUE,
                fixedColumns = TRUE,
                autoWidth = TRUE,
                ordering = TRUE,
                dom = 'tBp',
                buttons = list(
                  'excel', 'csv')
              ),
              class = "display",
              extensions = 'Buttons'
            ),
            server = TRUE
          )}

          output$anovaTable_flex_downloadData <- downloadHandler(
            filename = function() {
              paste("ANOVA_", Sys.Date(), ".xlsx", sep = "")
            },
            content = function(file) {
              # Write the full dataset to an Excel file
              writexl::write_xlsx(render_local_anova, file)
            }
          )

          render_local_tukey <-  if (tukey_needed) {
            local_tukey %>% tidy() %>%
              select(-null.value) %>%
              mutate(
                sig = case_when(
                  adj.p.value <= 0.001 ~ "***",
                  adj.p.value <= 0.01 ~ "**",
                  adj.p.value <= 0.05 ~ "*",
                  adj.p.value <= 0.1 ~ ".",
                  TRUE ~ ""
                )) %>%
              mutate(across(where(is.numeric), ~ round(., digits = 2)))} else {
                tibble::tibble()
              }
          output$tukeyTable_flex <- {
            DT::renderDT(
            DT::datatable(
              render_local_tukey,
              options = list(
                server = TRUE,
                filter = "top",
                selection = "multiple",
                searching = TRUE,
                fixedColumns = TRUE,
                autoWidth = TRUE,
                ordering = TRUE,
                dom = 'tBp',
                buttons = list(
                  'excel', 'csv')
              ),
              class = "display",
              extensions = 'Buttons'
            ),
            server = TRUE
          )
          }
          output$tukeyTable_flex_downloadData <- downloadHandler(
            filename = function() {
              paste("Tukey_", Sys.Date(), ".xlsx", sep = "")
            },
            content = function(file) {
              # Write the full dataset to an Excel file
              writexl::write_xlsx(render_local_tukey, file)
            }
          )



          render_local_letters <-  if (letters_needed) {
              local_names %>%
                rename(!!(
                  paste(combined_inputs$tukey_group, collapse = ":")
                ) := "group")
            } else {
                tibble::tibble()
            }
          output$tukeyLetters_flex  <- {
            DT::renderDT(
              DT::datatable(
                render_local_letters,
                options = list(
                  server = TRUE,
                  filter = "top",
                  selection = "multiple",
                  searching = TRUE,
                  fixedColumns = TRUE,
                  autoWidth = TRUE,
                  ordering = TRUE,
                  dom = 'tBp',
                  buttons = list(
                    'excel', 'csv')
                ),
                class = "display",
                extensions = 'Buttons'
              ),
              server = TRUE
            )
          }
          output$tukeyLetters_flex_downloadData <- downloadHandler(
            filename = function() {
              paste("Letters_", Sys.Date(), ".xlsx", sep = "")
            },
            content = function(file) {
              # Write the full dataset to an Excel file
              writexl::write_xlsx(render_local_letters, file)
            }
          )



        if (!combined_inputs$timeseries_plot) {
          set.seed(42)

          plot_d <- {
            local_table %>%
              ggplot(aes(
                x = as.POSIXct(timestamp_group, tz = 'UTC'),
                y = !!sym(combined_inputs$out_variables),
                color = letter,
                group = letter,
                fill = letter
              )) +
              geom_violin(
                na.rm = TRUE,
                width = 0.1,
                color = "black",
                alpha = 1,
                draw_quantiles = c(0.25, 0.5, 0.75),
                position = position_dodge(width = 0.8)
              ) +
              {
                if (nrow(local_table) > 2000) {
                  geom_jitter(
                    data = local_table %>%
                      {
                        if (all(most_interactive %in% colnames(.))){
                          . %>% group_by_at(most_interactive) %>%
                            slice_sample(n = 5) %>%
                            ungroup()
                        } else {
                          .
                        }
                      },
                    aes(
                      x = as.POSIXct(timestamp_group, tz = 'UTC'),
                      y = !!sym(combined_inputs$out_variables)
                    ),
                    height = 0,
                    width = 0.1,
                    color = "black",
                    alpha = 0.3
                  )
                } else {
                  geom_jitter(
                    height = 0,
                    width = 0.1,
                    color = "black",
                    alpha = 0.3
                  )
                }
              } +
              geom_label(
                aes(
                  label = letter,
                  y = quantile(
                    !!sym(combined_inputs$out_variables),
                    probs = seq(0, 1, .05),
                    na.rm = TRUE
                  )['100%'] * 1.01
                ),
                color = 'black',
                fill = 'white',
                alpha = 0.5,
                size = 4,
                vjust = 1,
                position = position_dodge(width = 0.8)
              ) +
              labs(x = 'time', y =
                     combined_inputs$out_variables) +
              scale_x_datetime() +
              scale_fill_viridis(discrete = TRUE) +
              facet_grid(combined_inputs$facet_formula,
                         labeller = label_both,
                         scales = 'free_x') +
              theme_gray() +
              theme(
                text = element_text(size = 12),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
              )
          }
        } else {
          set.seed(42)

          plot_d <- {
            local_table %>%
              ggplot(aes(
                x = as.POSIXct(timestamp_group, tz = 'UTC'),
                y = !!sym(combined_inputs$out_variables),
                color = !!sym(combined_inputs$gene_grouping),
                group = !!sym(combined_inputs$gene_grouping),
                fill = !!sym(combined_inputs$gene_grouping)
              )) +
              {
                if (nrow(local_table) > 2000) {
                  list(
                    geom_point(
                      data = local_table %>%
                        {
                          if (all(most_interactive %in% colnames(.))){
                            . %>% group_by_at(most_interactive) %>%
                              slice_sample(n = 5) %>%
                              ungroup()
                          } else {
                            .
                          }
                        },
                      aes(
                        x = as.POSIXct(timestamp_group, tz = 'UTC'),
                        y = !!sym(combined_inputs$out_variables)
                      ),
                      na.rm = TRUE,
                      width = 0.1,
                      alpha = 1,
                      draw_quantiles = c(0.25, 0.5, 0.75),
                      position = position_dodge(width = 0.8)
                    ),
                    geom_smooth(
                      data = local_table %>%
                        {
                          if (all(most_interactive %in% colnames(.))){
                            . %>% group_by_at(most_interactive) %>%
                              slice_sample(n = 5) %>%
                              ungroup()
                          } else {
                            .
                          }
                        },
                      aes(
                        x = as.POSIXct(timestamp_group, tz = 'UTC'),
                        y = !!sym(combined_inputs$out_variables)
                      )
                    )
                  )
                } else {
                  list(
                    geom_point(
                      na.rm = TRUE,
                      width = 0.1,
                      alpha = 1,
                      draw_quantiles = c(0.25, 0.5, 0.75),
                      position = position_dodge(width = 0.8)
                    ),
                    geom_smooth()
                  )
                }
              } +
              geom_label(
                aes(
                  label = letter,
                  y = quantile(
                    !!sym(combined_inputs$out_variables),
                    probs = seq(0, 1, .05),
                    na.rm = TRUE
                  )['100%'] * 1.01
                ),
                color = 'black',
                fill = 'white',
                alpha = 0.2,
                size = 4,
                vjust = 1,
                position = position_dodge(width = 0.8)
              ) +
              labs(x = 'time', y =
                     combined_inputs$out_variables) +
              scale_x_datetime(date_minor_breaks = "3 days") +
              scale_fill_viridis(discrete = TRUE) +
              facet_grid(combined_inputs$facet_formula,
                         labeller = label_both,
                         scales = 'free_x') +
              theme_gray() +
              theme(text = element_text(size = 12))
          }
        }

        output$distPlot <- renderPlot(execOnResize = FALSE, {
          plot_d
        })

        output$savePlot <- downloadHandler(
          filename = function() {
            "plot.svg"
          },
          content = function(file) {
            plot <- plot_d
            local_model <- local_model_d
            ggsave(
              file,
              plot + ggtitle(local_model),
              device = "svg",
              width = 16,
              height = 9
            )
          }
        )
      }#anova-tukey-letters-output

    },
    ignoreInit = FALSE,
    priority = 10
  )
}

# Run the application
shinyApp(ui = ui, server = server)
