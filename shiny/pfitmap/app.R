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

# We have problematic organisms, where multiple sequences of the same kind are assigned
# to the same taxon, probably always a species. Trying to get rid of those by filtering
# non-strain taxa from species with at least one strain.

# Step 1. Get all unique taxa
taxa = classified_proteins %>% 
  select(ncbi_taxon_id, tdomain:tstrain) %>% distinct()  %>% 
  mutate(tspecies = ifelse(is.na(tspecies) & ! is.na(tgenus), sprintf("%s sp.", tgenus), tspecies))

# Step 2. Left join with a list of species that have strains, and then filter.
taxa = taxa %>% 
  left_join(
    taxa %>% 
      filter(tspecies != tstrain) %>% 
      select(tdomain:tspecies) %>% distinct() %>% 
      mutate(strains=T),
    by = c("tdomain", "tkingdom", "tphylum", "tclass", "torder", "tfamily", "tgenus", "tspecies")
  ) %>% 
  replace_na(list('strains'=F)) %>%
  filter( ! ( strains & tspecies == tstrain ) )

# We will need a vector of protein superfamilies
psuperfamilies = (classified_proteins %>% select(psuperfamily) %>% distinct() %>% arrange(psuperfamily))$psuperfamily

# We also need a vector of databases
dbs = (classified_proteins %>% select(db) %>% distinct() %>% arrange(db))$db

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  sidebarLayout(
    sidebarPanel( 
      selectInput(
        'db', 'Database',
        dbs, selected = 'ref'
      ),
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
        ),
        selected = c('pclass')
      ),
      selectInput(
        'psuperfamilies', 'Protein superfamily',
        c('', psuperfamilies), selected = c('NrdGRE'), 
        multiple=T
      )
    ),
    mainPanel(
      h1('pfitmap'),
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
      filter(db == input$db) %>%
      inner_join(taxa %>% select(ncbi_taxon_id), by='ncbi_taxon_id') %>%
      filter(psuperfamily %in% dataPSuperfamiliesFilter()) %>%
      group_by_(input$taxonrank, input$proteinrank) %>%
      summarise(n=n())  %>% 
      spread_(input$proteinrank, 'n', fill=0)
  )
}

# Run the application 
shinyApp(ui = ui, server = server)