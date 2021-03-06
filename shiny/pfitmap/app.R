#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(Hmisc)       # Place early; it masks summarize from dplyr!
library(shiny)
library(readr)
library(dplyr)
library(data.table)
library(dtplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(stringr)
library(DT)
library(chorddiag)

# Some constants
PROTEIN_HIERARCHY = c( 'psuperfamily', 'pfamily', 'pclass', 'psubclass', 'pgroup' )
TAXON_HIERARCHY = c( 'tdomain', 'tkingdom', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies', 'tstrain' )

INDPROTEINS = 'indproteins'
COMBPROTEINS = 'combproteins'

DIV_PALETTE = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
DIV_PALETTE_96X = c(
  DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE
)
DIV_PALETTE_768X = c(
  DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X
)
DIV_PALETTE_6144X = c(
  DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X
)

LIGHT_PALETTE = c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99')
LIGHT_PALETTE_96X = c(
  LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE,
  LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE
)
LIGHT_PALETTE_768X = c(
  LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X
)
LIGHT_PALETTE_6144X = c(
  LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X
)
LIGHT_PALETTE_98304X = c(
  LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X,
  LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X
)

# Reading data and transforming
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

# We have problematic organisms, where multiple sequences of the same kind are
# assigned to the same taxon, a species or a genus. Trying to get rid of the
# species cases by filtering non-strain taxa from species with at least one
# strain. At the same time, delete all taxa with tgenus == tstrain and no
# tspecies.

# Step 1. Get all unique taxa
taxa = classified_proteins %>% 
  select(ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tstrain) %>% 
  distinct() %>% 
  filter( ! ( tgenus == tstrain & is.na(tspecies) ) ) %>%
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

# Fill in empty levels of the taxon hierarchy (can't be done before the steps
# involving taxa above).
classified_proteins = classified_proteins %>%
  mutate(
    tkingdom = ifelse(is.na(tkingdom), sprintf("%s, no kingdom", tdomain), tkingdom),
    tphylum = ifelse(is.na(tphylum), sprintf("%s, no phylum", sub(', no kingdom', '', tkingdom)), tphylum),
    tclass = ifelse(is.na(tclass), sprintf("%s, no class", sub(', no phylum', '', tphylum)), tclass),
    torder = ifelse(is.na(torder), sprintf("%s, no order", sub(', no class', '', tclass)), torder),
    tfamily = ifelse(is.na(tfamily), sprintf("%s, no family", sub(', no order', '', torder)), tfamily),
    tgenus = ifelse(is.na(tgenus), sprintf("%s, no genus", sub(', no family', '', tfamily)), tgenus),
    tspecies = ifelse(is.na(tspecies), sprintf("%s, no species", sub(', no genus', '', tgenus)), tspecies)
  )

# Do the same for taxa
taxa = taxa %>%
  mutate(
    tkingdom = ifelse(is.na(tkingdom), sprintf("%s, no kingdom", tdomain), tkingdom),
    tphylum = ifelse(is.na(tphylum), sprintf("%s, no phylum", sub(', no kingdom', '', tkingdom)), tphylum),
    tclass = ifelse(is.na(tclass), sprintf("%s, no class", sub(', no phylum', '', tphylum)), tclass),
    torder = ifelse(is.na(torder), sprintf("%s, no order", sub(', no class', '', tclass)), torder),
    tfamily = ifelse(is.na(tfamily), sprintf("%s, no family", sub(', no order', '', torder)), tfamily),
    tgenus = ifelse(is.na(tgenus), sprintf("%s, no genus", sub(', no family', '', tfamily)), tgenus),
    tspecies = ifelse(is.na(tspecies), sprintf("%s, no species", sub(', no genus', '', tgenus)), tspecies)
  )

# We will need a vector of protein superfamilies
psuperfamilies = (classified_proteins %>% select(psuperfamily) %>% distinct() %>% arrange(psuperfamily))$psuperfamily

# We also need a vector of databases
dbs = (classified_proteins %>% select(db) %>% distinct() %>% arrange(db))$db

# Define UI for application that draws a histogram
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel( 
      radioButtons(
        'protstattype', 'Type of protein statistic',
        list(
          'Individual proteins' = INDPROTEINS,
          'Combinations of proteins' = COMBPROTEINS
        ),
        selected = INDPROTEINS
      ),
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
      textOutput('debug'),
      tabsetPanel(type= 'tabs', 
        tabPanel('table', 
          fluidRow(
            column(4, checkboxInput('taxonomysort', 'Taxonomic sort', value=T)),
            column(4, uiOutput('trank4colour'))
          ),
          dataTableOutput('mainmatrix')
        ),
        tabPanel('chord graph',
          chorddiagOutput('chordgraph')
        ),
        tabPanel('distributions',
          selectInput(
            'sinastat', 'Statistic',
            list(
              'HMM score' = 'score', 'Sequence length' = 'seqlen', 'Alignment length' = 'align_length'
            )
          ),
          plotOutput('distgraph'))
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

    # Construct a field for taxonomical sort and on for tooltip for taxonomy.
    # This is done in two steps: first a string is constructed, then the the
    # string is used in a mutate_ statement.
    n = which(TAXON_HIERARCHY==input$taxonrank)
    ts_string = paste(
      'sprintf(strrep("%-50s",', which(TAXON_HIERARCHY == input$taxonrank),  '), ',
      paste(TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)], collapse=", "), ')'
    )
    ttt_string = paste(
      'paste(', paste(TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)], collapse=", "), ', sep="; ")'
    )
    t %>%
      mutate_(
        'tsort' = ts_string, 
        'tcolour' = input$trank4colour,
        'taxon_tooltip' = ttt_string
      )
  })

  # Returns a filtered and summarised table after applying the group by
  # criteria called for by the user.
  indproteins_sums_table = reactive({
    ###write(sprintf("indproteins_sums_table, protstattype: %s", input$protstattype), stderr())
    d = filtered_table() %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', input$taxonrank, input$proteinrank) %>%
      summarise(n=n()) %>%
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
      )
    
    d
  })

  # Calculates all present combinations of proteins at the proteinrank selected,
  # groups by the combinations, combines with n_genomes from taxa and returns.
  combproteins_sums_table = reactive({
    ###write(sprintf("combproteins_sums_table, protstattype: %s", input$protstattype), stderr())
    d = filtered_table() %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', 'ncbi_taxon_id', input$taxonrank, input$proteinrank) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      spread_(input$proteinrank, input$proteinrank, fill='') %>%
      select(-n)
    d = d %>%
      unite(comb, 6:length(colnames(d)), sep=':') %>%
      mutate(comb=sub('::*', ':', sub(':*$', '', sub('^:*', '', comb)))) %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', input$taxonrank, 'comb') %>%
      summarise(n=n()) %>%
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
      )

    d
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

  output$trank4colour = renderUI({
    ranks = list()
    for ( r in TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)] ) { 
      ranks[[capitalize(sub('^t', '', r))]] = r 
    }
    selectInput(
      'trank4colour', 'Colour by taxon',
      ranks, selected = 'tdomain'
    )
  })
  
  output$mainmatrix = renderDataTable(
    {
      t = switch(
        input$protstattype,
        indproteins  = indproteins_sums_table() %>% spread_(input$proteinrank, 'n', fill=0),
        combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0)
      )
      ###write(sprintf("colnames(t): %s", colnames(t)), stderr())
      
      if ( input$taxonomysort ) {
        t = t %>% arrange(tsort)
      }
      
      # This is to get the right column names, a bit involved perhaps...
      t = t %>% mutate_('Taxon'=input$taxonrank, `N. genomes`='n_genomes') %>%
        mutate(Taxon = sprintf("<span title='%s'>%s</span>", taxon_tooltip, Taxon))
      c = colnames(t)
      ###write(sprintf("c: %s", c), stderr())
      t = t %>%
        select(tcolour, c(length(c)-1,length(c),8:length(c)-2))
      datatable(
        t, 
        rownames=F, 
        escape = c(T, F),
        options=list(
          lengthMenu = c(50, 100, 250, 500),
          columnDefs = list(list(targets = 0, visible = FALSE))
        )
      ) %>%
        formatStyle(
          'Taxon', 'tcolour',
          backgroundColor = styleEqual(
            unique(t$tcolour), LIGHT_PALETTE_98304X[1:length(unique(t$tcolour))]
          )
        )
    }
  )
  
  output$chordgraph = renderChorddiag({
    t = switch(
      input$protstattype,
      indproteins  = indproteins_sums_table() %>% spread_(input$proteinrank, 'n', fill=0),
      combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0)
    ) %>% 
      select(-tsort, -tcolour, -taxon_tooltip, -n_genomes)
    write(sprintf("colnames(t): %s", colnames(t)), stderr())
    m = as.matrix(t[,2:length(colnames(t))])
    rownames(m) = (t %>% select(t=1))$t
    chorddiag(m, type = "bipartite")
  })
  
  output$distgraph = renderPlot({
    subc =  ifelse(
      which(PROTEIN_HIERARCHY==input$proteinrank) == length(PROTEIN_HIERARCHY),
      input$proteinrank, PROTEIN_HIERARCHY[which(PROTEIN_HIERARCHY==input$proteinrank) + 1]
    )
    d = filtered_table() %>%
      mutate(seqlen = str_length(seq)) %>%
      mutate_(
        'stat' = input$sinastat,
        'c' = input$proteinrank,
        'subc' = subc
      )
    
    ggplot(d, aes(x=c, y=stat)) + 
      geom_violin() +
      geom_sina(aes(colour=subc), method='counts') +
      scale_colour_manual(sprintf('Protein %s', sub('^p', '', subc)), values=DIV_PALETTE_768X) +
      xlab(sprintf("Protein %s", sub('^p', '', input$proteinrank))) +
      ylab('Statistic')
  }) 
  
  output$ssversion = renderText(
    (classified_proteins %>% 
      transmute(ssversion = sprintf("Source database: %s %s, downloaded %s", ss_source, ss_name, ss_version)) %>% 
      distinct())$ssversion
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
