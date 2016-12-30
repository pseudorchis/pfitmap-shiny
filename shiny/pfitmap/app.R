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
library(data.table)
library(dtplyr)
library(tidyr)

classified_proteins = data.table(
  read_tsv(
    Sys.getenv('PFITMAP_DATA'),
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
)

# We will need a vector of protein superfamilies
psuperfamilies = (classified_proteins %>% select(psuperfamily) %>% distinct() %>% arrange(psuperfamily))$psuperfamily

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  sidebarLayout(
    sidebarPanel( 
      selectInput(
        'taxonrank', 'Taxon rank', 
        list(
          'Domain'  = 'tdomain',  'Phylum'  = 'tphylum', 'Class'   = 'tclass',
          'Order'   = 'torder',   'Family'  = 'tfamily', 'Genus'   = 'tgenus',
          'Species' = 'tspecies', 'Strain'  = 'tstrain'
        )
      ),
      selectInput(
        'proteinrank', 'Protein rank',
        list(
          'Superfamily' = 'psuperfamily', 'Family'   = 'pfamily',
          'Class'       = 'pclass',       'Subclass' = 'psubclass',
          'Group'       = 'pgroup'
        )
      ),
      selectInput(
        'psuperfamilies', 'Protein superfamily',
        c('', psuperfamilies), multiple=T
      )
    ),
    mainPanel(
      h1('pfitmap'),
      textOutput('psupdebug'),
      dataTableOutput('mainmatrix')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  dataPSuperfamiliesFilter <- reactive({
    if ( length(input$psuperfamilies) == 0 ) {
      return(psuperfamilies)
    }
    return(input$psuperfamilies)
  })
  
  output$mainmatrix = renderDataTable(
    classified_proteins %>% 
      filter(psuperfamily %in% dataPSuperfamiliesFilter()) %>%
      group_by_(input$taxonrank, input$proteinrank) %>%
      summarise(n=n())  %>% 
      spread_(input$proteinrank, 'n', fill=0)
  )
}

# Run the application 
shinyApp(ui = ui, server = server)