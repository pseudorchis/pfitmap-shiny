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
library(ggplot2)
library(ggforce)

TAXON_HIERARCHY = c( 'tdomain', 'tkingdom', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies', 'tstrain' )

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
      score = col_double()
    )
  ) %>%
    mutate(
      pfamily = ifelse(is.na(pfamily), sprintf("%s, no family", psuperfamily), pfamily),
      pclass = ifelse(is.na(pclass), sprintf("%s, no class", pfamily), pclass),
      psubclass = ifelse(is.na(psubclass), sprintf("%s, no subclass", pclass), psubclass),
      pgroup = ifelse(is.na(pgroup), sprintf("%s, no group", psubclass), pgroup)
    )
)

# We have problematic organisms, where multiple sequences of the same kind are assigned
# to the same taxon, probably always a species. Trying to get rid of those by filtering
# non-strain taxa from species with at least one strain.

# Step 1. Get all unique taxa
taxa = classified_proteins %>% 
  select(ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tstrain) %>% 
  distinct() %>% 
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
      wellPanel(
        selectInput(
          'psuperfamilies', 'Protein superfamilies',
          c('', psuperfamilies), selected = c('NrdGRE'), 
          multiple=T
        ),
        uiOutput('pfamilies'),
        uiOutput('pclasses')
      )
    ),
    mainPanel(
      h1('pfitmap'),
      textOutput('ssversion'),
      tabsetPanel(type= 'tabs', 
        tabPanel('table', 
          checkboxInput('taxonomysort', 'Taxonomic sort'),
          dataTableOutput('mainmatrix')
        ),
        tabPanel('graph',  plotOutput('maingraph'))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Reactive methods outside assignments
  
  # Returns a table after applying all filters the user have called for
  filtered_table = reactive({
    t = classified_proteins %>% 
      filter(db == input$db) %>%
      inner_join(taxa %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    
    if ( length(input$psuperfamilies) > 0 ) {
      t = t %>% filter(psuperfamily %in% input$psuperfamilies)
    }
    
    if ( length(input$pfamilies) > 0 ) {
      t = t %>% filter(pfamily %in% input$pfamilies)
    }
    
    if ( length(input$pclasses) > 0 ) {
      t = t %>% filter(pclass %in% input$pclasses)
    }

    t
  })

  # Returns a filtered and summarised table after applying the group by
  # criteria called for by the user.
  group_summed_table = reactive({
    filtered_table() %>%
      group_by_(input$taxonrank, input$proteinrank) %>%
      summarise(n=n())  %>% 
      inner_join(
        taxa %>%
          inner_join(
            classified_proteins %>%
              filter(db == input$db) %>%
              select(ncbi_taxon_id) %>% distinct(),
            by = 'ncbi_taxon_id'
          ) %>%
          group_by_(input$taxonrank) %>%
          summarise(n_genomes = n()),
        by=c(input$taxonrank)
      ) %>%
      spread_(input$proteinrank, 'n', fill=0)
  })

  output$pfamilies = renderUI({
    pf = classified_proteins
    if ( length(input$psuperfamilies) > 0 ) {
      pf = pf %>% filter(psuperfamily %in% input$psuperfamilies)
    }
    pf = (pf %>% select(pfamily) %>% distinct())$pfamily
    
    selectInput(
      'pfamilies', 'Protein families',
      pf, multiple = T
    )
  })
  
  output$pclasses = renderUI({
    pc = classified_proteins
    if ( length(input$pfamilies) > 0) {
      pc = pc %>% filter(pfamily %in% input$pfamilies)
    }
    pc = (pc %>% select(pclass) %>% distinct() %>% arrange(pclass))$pclass
    
    selectInput(
      'pclasses', 'Protein classes',
      pc, multiple = T
    )
  })
  
  output$mainmatrix = renderDataTable({
    t = group_summed_table()
    # This is to get the right column names, a bit involved perhaps...
    c = colnames(t)
    t %>% mutate_('Taxon'=input$taxonrank, `N. genomes`='n_genomes') %>% 
      select(c(length(c)+1,length(c)+2,3:length(c)))
  })
  
  output$maingraph = renderPlot({
    ggplot(filtered_table(), aes(x=pclass, y=score, color=psubclass)) + geom_sina()
  }) 
  
  output$ssversion = renderText(
    (classified_proteins %>% 
      transmute(ssversion = sprintf("Source database: %s %s, downloaded %s", ss_source, ss_name, ss_version)) %>% 
      distinct())$ssversion
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
