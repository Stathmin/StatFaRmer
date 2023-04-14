library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(here)
library(factoextra)

`%>%` <- magrittr::`%>%`
set.seed(42)

# Load dataset
datatable <-
  readRDS(paste0(here::here(), '/.cache/printable_table.rds')) %>%
  tidyr::pivot_wider(
    names_from = c(tukey_grouping, grouping_gene),
    values_from = c(tukey_group, tukey_letter),
    names_sep = '.'
  ) %>%
  dplyr::select(-timestamp) %>%
  dplyr::distinct() %>%
  dplyr::rename_with(\(x) stringi::stri_replace_all_regex(x,
                                                          'tukey_group\\.',
                                                          '')) %>%
  dplyr::select(-dplyr::starts_with('treatment.')[-1]) %>%
  dplyr::select(-dplyr::starts_with('tukey_letter.treatment.')[-1]) %>%
  dplyr::rename_with(\(x) {
    stringi::stri_replace_all_regex(x,
                                    '(?<=^(treatment\\.|tukey_letter\\.treatment\\.)).*',
                                    '')
  }) %>%
  dplyr::mutate(., dplyr::across(-c(trait_value), as.character)) %>%
  tidyr::pivot_wider(
    names_from = trait,
    values_from = c(dplyr::contains('letter.'), trait_value),
    names_sep = '.'
  ) %>%
  dplyr::rename_with(\(x) {
    stringi::stri_replace_all_regex(x,
                                    '\\.{2,}',
                                    '\\.')
  }) %>%
  dplyr::select(dplyr::where(\(x) dplyr::n_distinct(x) > 1)) %>%
  dplyr::rename_with(\(x) stringi::stri_replace_all_fixed(x,
                                                          'trait_value.', '')) %>%
  dplyr::rename_with(\(x) stringi::stri_replace_all_fixed(x,
                                                          '_logit', '')) %>%
  dplyr::relocate(treatment., .before = var) %>%
  dplyr::select(-dplyr::matches('bin[0125]'))

datatable %>% write.csv('test.csv')

non_numeric_columns <- datatable %>%
  dplyr::select(!is.numeric) %>%
  dplyr::select(dplyr::where(\(x) dplyr::n_distinct(x) < length(x) / 2)) %>%
  colnames()

numeric_columns <- datatable %>%
  dplyr::select(is.numeric) %>%
  colnames()

res.pca <- datatable %>%
  dplyr::select(numeric_columns) %>%
  prcomp(scale. = TRUE)

# Define UI for app that draws a histogram ----
ui <-  function(req) {
  fluidPage(# App title ----
            titlePanel("Please Call an Ambulance"),

            # Sidebar layout with input and output definitions ----
            sidebarLayout(
              # Sidebar panel for inputs ----
              sidebarPanel(
                # Input: Selectors ----

                selectInput(
                  inputId = "grouping_column",
                  label = "Grouping column",
                  choices = non_numeric_columns
                ),

                selectInput(
                  inputId = "vectors",
                  label = "Vectors",
                  multiple = TRUE,
                  choices = numeric_columns
                )

              ),

              # Main panel for displaying outputs ----
              mainPanel(# Output: Histogram ----
                        plotOutput(outputId = "PCA"),)
            ))
}


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  groups <- reactive({
    datatable %>%
      dplyr::select(dplyr::one_of(input$grouping_column)) %>%
      dplyr::pull() %>%
      as.character()
  })

  vectors <- reactive({
    input$vectors
  })

  output$PCA <- renderPlot({
    if (length(vectors()) > 0) {
      fviz_pca_biplot(
        res.pca,
        habillage = groups(),
        # color by groups
        addEllipses = TRUE,
        # Concentration ellipses
        ellipse.type = "confidence",
        ellipse.level = 0.95,
        legend.title = "Groups",
        repel = TRUE,
        select.var = list(names = vectors()),
        label = "var"
      )
    } else {
      fviz_pca_ind(
        res.pca,
        habillage = groups(),
        # color by groups
        addEllipses = TRUE,
        # Concentration ellipses
        ellipse.type = "confidence",
        ellipse.level = 0.95,
        legend.title = "Groups",
        repel = TRUE,
        label = "none"
      )
    }
  })
}

shinyApp(ui = ui, server = server)
