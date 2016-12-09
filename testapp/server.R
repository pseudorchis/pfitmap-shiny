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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  cp = read_tsv(
    '../data/test/classified_proteins.10000.tsv',
    col_types = cols(
      .default = col_character(),
      profile_length = col_integer(), align_length = col_integer(),
      align_start = col_integer(), align_end = col_integer(),
      prop_matching = col_double(), e_value = col_double(),
      score = col_double()
    )
  )
  genomes = cp %>% select(tdomain:tstrain) %>% distinct()
  
  tdomain_pfamily <- genomes %>%
    group_by(tdomain) %>%
    summarise(ngenomes = n()) %>%
    ungroup() %>%
    inner_join(
      cp %>%
        group_by(tdomain, psuperfamily, pfamily) %>%
        summarise(nproteins=n()) %>%
        ungroup(),
      by=c('tdomain')
    )
  
  output$distPlot <- renderPlot({
    tphylum_pfamily %>% 
      transmute(
        taxon=tdomain,
        protein=sprintf("%s:%s", psuperfamily, pfamily), 
        ngenomes, nproteins
      ) %>% 
      ggplot(aes(x=taxon, fill=protein)) + geom_bar()
  })
})
