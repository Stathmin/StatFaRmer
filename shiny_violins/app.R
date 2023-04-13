library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(here)

`%>%` <- magrittr::`%>%`
set.seed(42)

thematic::thematic_shiny(font = "auto")

# Load dataset
  datatable <- readRDS(paste0(here::here(),'/.cache/printable_table.rds'))
  traits <- datatable %>%
    dplyr::distinct(trait)
  tukey_groupings <- datatable %>%
    dplyr::distinct(tukey_grouping)
  grouping_genes <- datatable %>%
    dplyr::distinct(grouping_gene)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  theme = bs_theme(
    bg = "#333", fg = "white", primary = "#FCC780",
    base_font = font_google("Fira Code"),
    code_font = font_google("Fira Code")
  ),

  # App title ----
  titlePanel("Violins Apocalyptica"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Selectors ----
      selectInput(inputId = "traits",
                  label = "Trait",
                  choices = traits),
      selectInput(inputId = "tukey_groupings",
                  label = "Factor",
                  choices = tukey_groupings),
      selectInput(inputId = "grouping_genes",
                  label = "Group of genes",
                  choices = grouping_genes)

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot"),

      #Output: Table ----
      DT::dataTableOutput('table')

    )
  )
)



# Define server logic required to draw a histogram ----
server <- function(input, output) {

  local_datatable <- reactive({datatable %>%
    dplyr::filter( trait          == input$traits,
                   tukey_grouping == input$tukey_groupings,
                   grouping_gene  == input$grouping_genes)
  })

  output$distPlot <- renderPlot({
    ggplot(
      local_datatable(),
      aes(y = trait_value,
                   x = tukey_group,
                   color = tukey_letter)
    ) +
      geom_violin(
        draw_quantiles = TRUE,
        scale = 'count'
      ) +
      geom_boxplot(width = 0.1)

  })
  output$table <- DT::renderDataTable(
    local_datatable() %>%
      dplyr::select(
        dplyr::where(\(x) dplyr::n_distinct(x) > 1)
        ) %>%
      dplyr::mutate(
        trait_value = round(trait_value, digits = 2)
      )
    ,
    options = list(lengthChange = FALSE)
  )

}

shinyApp(ui = ui, server = server)




