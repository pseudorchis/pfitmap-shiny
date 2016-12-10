#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

# Constants

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  dataMatrix <- reactive({
    data <- genomes %>%
      group_by_(.dots = c(as.symbol(input$taxon_rank))) %>%
      summarise(ngenomes = n()) %>%
      ungroup() %>%
      inner_join(
        cp %>%
          filter(psuperfamily %in% dataPSuperfamiliesFilter()) %>%
          group_by_(.dots = lapply(c(input$taxon_rank, input$protein_rank), as.symbol)) %>%
          summarise(nproteins=n()) %>%
          ungroup(),
        by=c(input$taxon_rank)
      ) %>%
      mutate_('taxon' = input$taxon_rank, 'protein' = input$protein_rank)
  })
  
  dataPSuperfamiliesFilter <- reactive({
    if ( length(input$psuperfamilies) == 0 ) {
      return(psuperfamilies$psuperfamily)
    }
    input$psuperfamilies
  })

  output$matrix <- renderTable({
    dataMatrix() %>%
      select(taxon, protein, ngenomes, nproteins) %>%
      spread(protein, nproteins)
  })
  
  output$summaryplot <- renderPlot({
    ggplot(dataMatrix(), aes(x=taxon, y=nproteins, colour=protein)) +
      geom_point() +
      theme(
        axis.text.x = element_text(angle=60, hjust=1)
      ) +
      scale_y_log10()
  })
})
