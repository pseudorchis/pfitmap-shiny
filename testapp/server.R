#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Constants

# Initialize tables from disk
cp <- read_tsv(
  '../data/test/classified_proteins.10000.tsv',
  col_types = cols(
    .default = col_character(),
    profile_length = col_integer(), align_length = col_integer(),
    align_start = col_integer(), align_end = col_integer(),
    prop_matching = col_double(), e_value = col_double(),
    score = col_double()
  )
) %>%
  mutate(
    tdomain = ifelse(is.na(tdomain), 'No domain', tdomain),
    tkingdom = ifelse(is.na(tkingdom), sprintf('%s unsp.', tdomain), tkingdom),
    tphylum = ifelse(is.na(tphylum), sprintf('%s unsp.', tdomain), tphylum),
    tclass = ifelse(is.na(tclass), sprintf('%s unsp.', tphylum), tclass),
    torder = ifelse(is.na(torder), sprintf('%s unsp.', tclass), torder),
    tfamily = ifelse(is.na(tfamily), sprintf('%s unsp.', torder), tfamily),
    tgenus = ifelse(is.na(tgenus), sprintf('%s unsp.', tfamily), tgenus),
    tspecies = ifelse(is.na(tspecies), sprintf('%s unsp.', tgenus), tspecies),
    pfamily = ifelse(is.na(pfamily), sprintf("%s unsp.", psuperfamily), pfamily),
    pclass = ifelse(is.na(pclass), sprintf("%s unsp.", pfamily), pclass),
    psubclass = ifelse(is.na(psubclass), sprintf("%s unsp.", pclass), psubclass),
    pgroup = ifelse(is.na(pgroup), sprintf("%s unsp.", psubclass), pgroup)
  )
genomes <- cp %>% select(tdomain:tstrain) %>% distinct()

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$distPlot <- renderPlot({
    genomes %>%
      group_by_(.dots = c(as.symbol(input$taxon_rank))) %>%
      summarise(ngenomes = n()) %>%
      ungroup() %>%
      inner_join(
        cp %>%
          group_by_(.dots = lapply(c(input$taxon_rank, input$protein_rank), as.symbol)) %>%
          summarise(nproteins=n()) %>%
          ungroup(),
        by=c(input$taxon_rank)
      ) %>%
      mutate_('taxon' = input$taxon_rank, 'protein' = input$protein_rank) %>%
      ggplot(aes(x=taxon, y=nproteins, fill=protein)) + 
        geom_bar(stat='identity') +
        theme(
          axis.text.x = element_text(angle=60, hjust=1)
        )
  })
})
