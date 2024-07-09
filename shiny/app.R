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

options(shiny.reactlog = TRUE)
set.seed(42)

formulate <- function(RHS, LHF) {
  if (length(LHF) == 0) {
    LHS = '1'
  } else if ((length(LHF) == 1)) {
    LHS = str_interp('1 + ${LHF}')
  } else {
    LHS = str_interp('1 + (${paste(LHF, collapse=" + ")})^3')
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
        ,
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
        selected = named_timestamp_groups[c(1,
                                            round(
                                              length(named_timestamp_groups) /2
                                              ),
                                            length(named_timestamp_groups))] %>%
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
        'Facet formula',
        value = 'treatment ~ timestamp_group',
        placeholder = 'treatment ~ timestamp_group'
      ),
      checkboxInput('deltas', ': use deltas', FALSE)

    ),

    mainPanel(
      uiOutput("formula"),
      plotOutput("distPlot", height = "800px", width = '1200px'),
      tabsetPanel(
        tabPanel(
          "Raw Table",
          uiOutput("RAW"),
          verbatimTextOutput("raw"),
          fluidRow(uiOutput("raw_flex"))
        ),
        tabPanel(
          "Descriptive",
          uiOutput("DESCvarNames"),
          verbatimTextOutput("DescriptiveTable"),
          fluidRow(uiOutput("descriptiveTable_flex"))
        ),
        tabPanel(
          "ANOVA",
          uiOutput("ANOVAvarNames"),
          verbatimTextOutput("anovaTable"),
          fluidRow(uiOutput("anovaTable_flex"))
        ),
        tabPanel(
          "Tukey",
          uiOutput("TUKEYvarNames"),
          verbatimTextOutput("tukeyTable"),
          fluidRow(uiOutput("tukeyTable_flex"))
        ),
        tabPanel(
          "Group Letters",
          uiOutput("TUKEYLvarNames"),
          verbatimTextOutput("tukeyLetters"),
          fluidRow(uiOutput("tukeyLetters_flex"))
        )
      )
    )
  ))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  set.seed(42)

  updateFlag <- reactiveVal(TRUE)

  gene_grouping_d <- reactive({
    input$gene_grouping
  })

  observeEvent(gene_grouping_d(), {

    updateFlag(FALSE)

    isolate({
      updateSelectizeInput(
        session,
        'gene_groups',
        choices = c(),
        selected = c(),
        server = TRUE
      )
    })

    unique_groupings <- get_unique(merged_table, gene_grouping_d())
    updateSelectizeInput(
      session,
      'gene_groups',
      choices = unique_groupings,
      selected = unique_groupings,
      server = TRUE
    )

    updateFlag(TRUE)

  }, ignoreInit = TRUE)

  gene_groups_d <- reactive({
    req(isTRUE(updateFlag()))
    {
      input$gene_groups
    }
  })
  treatments_d <- reactive(input$treatments) %>% debounce(1000)
  cultivars_d <- reactive(input$cultivars) %>% debounce(1000)
  timestamp_groups_d <-
    reactive(unname(input$timestamp_groups)) %>% debounce(1000)
  out_variables_d <- reactive(input$out_variables)

  grouping_factors_d <- reactive({
    sort(input$grouping_factors)
  }) %>% debounce(1000)

  local_model_d <- reactive({
    formulate(out_variables_d(), grouping_factors_d())
  })
  tukey_group_d <- reactive({
    sort(input$tukey_group)
  })

  initial_table <- reactive({
    req(gene_groups_d())
    req(gene_grouping_d())
    req(isTRUE(updateFlag()))

    merged_table %>%
      filter(
        treatment %in% treatments_d(),
        cultivar %in% cultivars_d(),
        as.character(timestamp_group) %in% as.character(timestamp_groups_d()),
        !!sym(gene_grouping_d()) %in% gene_groups_d()
      ) %>%
      arrange(timestamp_group) %>%
      mutate(timestamp_group = as.factor(timestamp_group)) %>%
      unite(group,
            !!sort(tukey_group_d()),
            sep = ":",
            remove = FALSE) %>%
      select(-c(where(is.numeric), -!!out_variables_d()))
  })
  current_table <- reactive(if (input$deltas) {
    req(isTRUE(updateFlag()))

    initial_table() %>%
      arrange(dbscan_cluster) %>%
      group_by(v_t_r, dbscan_cluster) %>%
      mutate(across(any_of(out_variables_d()), \(x) (median(x, na.rm = TRUE)))) %>%
      ungroup() %>%
      select(-timestamp) %>%
      distinct() %>%
      group_by(v_t_r) %>%
      mutate(across(any_of(out_variables_d()), \(x) (x - lag(x)))) %>%
      ungroup() %>%
      group_by(cultivar) %>%
      mutate(across(any_of(out_variables_d()), \(x) (x / x[treatment == 1]))) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), ~ na_if(., Inf)), across(where(is.numeric), ~
                                                                  na_if(., -Inf))) %>%
      drop_na(out_variables_d())
  } else {
    initial_table()
  })

  {
    local_descriptive <- reactive({
      current_table() %>%
        filter(!is.na(!!sym(out_variables_d()))) %>%
        arrange(!!sym(out_variables_d())) %>%
        unite(group,
              !!tukey_group_d(),
              sep = ":",
              remove = FALSE) %>%
        select(group, !!out_variables_d()) %>%
        group_by(group) %>%
        summarise(
          n = dplyr::n(),
          median = median(!!rlang::sym(out_variables_d())),
          mean = mean(!!rlang::sym(out_variables_d())),
          cv_perc = 100 * sd(!!rlang::sym(out_variables_d())) / median(!!rlang::sym(out_variables_d())),
          min = min(!!rlang::sym(out_variables_d())),
          max = max(!!rlang::sym(out_variables_d())),
          skewness = moments::skewness(!!rlang::sym(out_variables_d())),
          kurtosis = moments::kurtosis(!!rlang::sym(out_variables_d())),
        ) %>%
        ungroup() %>%
        arrange(mean)
    })
    local_anova <- reactive({
      aov(as.formula(local_model_d()), data = current_table())
    })
    tukey_needed <- reactive({
      as.formula(local_model_d()) %>%
        all.vars() %>%
        length() > 1
    })
    letters_needed <- reactive({
      tukey_needed() &
        length(vecsets::vsetdiff(sort(tukey_group_d()), sort(grouping_factors_d()))) == 0
    })
    local_tukey <- reactive(if (tukey_needed()) {
      TukeyHSD(local_anova(), ordered = TRUE) %>%
        purrr::modify_depth(5, \(x) replace_na(x, replace = 0), .ragged = TRUE)
    } else {
      NA
    })
    local_names <- reactive(if (letters_needed()) {
      multcompLetters4(local_anova(), local_tukey()) %>%
        .[[paste(sort(tukey_group_d()),
                 sep = ':',
                 collapse = ':')]] %>%
        as.data.frame.list() %>%
        as_tibble(rownames = 'group') %>%
        select('group', 'Letters') %>%
        rename(letter = Letters)
    } else {
      NA
    })
    local_table <- reactive(if (letters_needed()) {
      current_table() %>%
        filter(!is.na(!!out_variables_d())) %>%
        left_join(local_names(), by = c('group' = 'group'))
    } else {
      current_table() %>%
        mutate(letter = '-')
    })
  }#anova-tukey-letters-calc

  {
    output$formula <- renderUI({
      HTML(as.character(
        div(style = "text-align: center; font-weight: bold;", local_model_d())
      ))
    })
    output$descriptiveTable_flex <- renderUI({
      local_descriptive()
      local_descriptive()  %>%
        mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>%
        DT::datatable(
          filter = "top",
          selection = "multiple",
          escape = FALSE,
          style = 'auto',
          extensions = 'Buttons',
          options = list(
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            ordering = TRUE,
            dom = 'tBp',
            buttons = list((
              list(
                extend = "excel",
                text = "Download Full Results",
                filename = 'descriptive',
                exportOptions = list(modifier = list(page = "all"))
              )
            ))
          ),
          class = "display"
        ) %>%
        DT::renderDataTable(server = FALSE) %>%
        fluidPage()
    })
    output$anovaTable_flex <- renderUI({
      local_anova()
      tidy(local_anova()) %>%
        mutate(
          sig = case_when(
            p.value <= 0.001 ~ "***",
            p.value <= 0.01 ~ "**",
            p.value <= 0.05 ~ "*",
            p.value <= 0.1 ~ ".",
            .default = ""
          )
        )  %>%
        mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>%
        DT::datatable(
          filter = "top",
          selection = "multiple",
          escape = FALSE,
          style = 'auto',
          extensions = 'Buttons',
          options = list(
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            ordering = TRUE,
            dom = 'tBp',
            buttons = list((
              list(
                extend = "excel",
                text = "Download Full Results",
                filename = 'anova',
                exportOptions = list(modifier = list(page = "all"))
              )
            ))
          ),
          class = "display"
        ) %>%
        DT::renderDataTable(server = FALSE) %>%
        fluidPage()
    })
    output$tukeyTable_flex <- renderUI(if (tukey_needed()) {
      tidy(local_tukey()) %>%
        select(-null.value) %>%
        mutate(
          sig = case_when(
            adj.p.value <= 0.001 ~ "***",
            adj.p.value <= 0.01 ~ "**",
            adj.p.value <= 0.05 ~ "*",
            adj.p.value <= 0.1 ~ ".",
            .default = ""
          ),
          across(where(is.numeric), \(x) round(x, digits = 2))
        ) %>%
        DT::datatable(
          filter = "top",
          selection = "multiple",
          escape = FALSE,
          style = 'auto',
          extensions = 'Buttons',
          options = list(
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            ordering = TRUE,
            dom = 'tBp',
            buttons = list((
              list(
                extend = "excel",
                text = "Download Full Results",
                filename = 'tukey',
                exportOptions = list(modifier = list(page = "all"))
              )
            ))
          ),
          class = "display"
        ) %>%
        DT::renderDataTable(server = FALSE) %>%
        fluidPage()
    } else {
      ''
    })
    output$tukeyLetters_flex <- renderUI(if (letters_needed()) {
      {
        local_names() %>%
          rename(!!(paste(tukey_group_d(), collapse = ":")) := "group")
      } %>%
        mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>%
        DT::datatable(
          filter = "top",
          selection = "multiple",
          escape = FALSE,
          style = 'auto',
          extensions = 'Buttons',
          options = list(
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            ordering = TRUE,
            dom = 'tBp',
            buttons = list((
              list(
                extend = "excel",
                text = "Download Full Results",
                filename = 'letters',
                exportOptions = list(modifier = list(page = "all"))
              )
            ))
          ),
          class = "display"
        ) %>%
        DT::renderDataTable(server = FALSE) %>%
        fluidPage()
    } else {
      ''
    })
    output$raw_flex <- renderUI({
      local_table()
      local_table() %>%
        mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>%
        DT::datatable(
          filter = "top",
          selection = "multiple",
          escape = FALSE,
          style = 'auto',
          extensions = 'Buttons',
          options = list(
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            ordering = TRUE,
            dom = 'tBp',
            buttons = list((
              list(
                extend = "excel",
                text = "Download Full Results",
                filename = 'raw',
                exportOptions = list(modifier = list(page = "all"))
              )
            ))
          ),
          class = "display"
        ) %>%
        DT::renderDataTable(server = FALSE) %>%
        fluidPage()
    })
  }#anova-tukey-letters-output

  facet_formula_d <-
    reactive(as.formula(input$facet_formula)) %>% debounce(3000)

  output$distPlot <- renderPlot(execOnResize = FALSE, {
    local_table() %>%
      ggplot(aes(
        x = as.POSIXct(timestamp_group, tz = 'UTC'),
        y = !!sym(out_variables_d()),
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
      geom_jitter(
        height = 0,
        width = 0.1,
        color = "black",
        alpha = 0.3
      ) +
      geom_label(
        aes(
          label = letter,
          y = quantile(
            !!sym(out_variables_d()),
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
      labs(x = 'time',
           y = ifelse(
             input$deltas,
             paste0('delta_', out_variables_d()),
             out_variables_d()
           )) +
      scale_x_datetime() +
      scale_fill_viridis(discrete = TRUE) +
      facet_grid(facet_formula_d(),
                 labeller = label_both,
                 scales = 'free_x') +
      theme_gray() +
      theme(
        text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

  })
}

# Run the application
shinyApp(ui = ui, server = server)
