#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(tidyr)

classified_proteins = read_tsv(
  '../../../data/test/classified_proteins.10000.tsv',
  col_types = cols(
    .default = col_character(),
    profile_length = col_integer(),
    align_length = col_integer(),
    align_start = col_integer(),
    align_end = col_integer(),
    prop_matching = col_double(),
    ss_version = col_integer(),
    e_value = col_double(),
    score = col_double(),
    profile_version = col_double()
  )
)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Eriks och Daniels testapp"),
   
  fluidRow(
    column(2),
    column(4, 
      selectInput(
        'taxonrank', 'Taxon rank', 
        list(
          'Domain'  = 'tdomain',  'Phylum'  = 'tphylum', 'Class'   = 'tclass',
          'Order'   = 'torder',   'Family'  = 'tfamily', 'Genus'   = 'tgenus',
          'Species' = 'tspecies', 'Strain'  = 'tstrain'
        )
      )
    ),
    column(4,
      selectInput(
        'proteinrank', 'Protein rank',
        list(
          'Superfamily' = 'psuperfamily', 'Family'   = 'pfamily',
          'Class'       = 'pclass',       'Subclass' = 'psubclass',
          'Group'       = 'pgroup'
        )
      )
    )
  ),
  
  fluidRow(
    column(12, dataTableOutput('mainmatrix'))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$mainmatrix = renderDataTable(
    classified_proteins %>% group_by_(input$taxonrank, input$proteinrank) %>%
      summarise(n=n())  %>% 
      spread_(input$proteinrank, 'n', fill=0)
  )
}

# Run the application 
shinyApp(ui = ui, server = server)