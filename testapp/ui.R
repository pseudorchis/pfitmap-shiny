#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("Pfitmap testapp"),
    
    # Dropboxes for grouping at the top
    fluidRow(
      column(
        6,
        selectInput(
          'taxon_rank', 'Taxon rank', 
          c(
            'Domain' = 'tdomain', 'Kingdom' = 'tkingdom', 'Phylum' = 'tphylum',
            'Class' = 'tclass', 'Order' = 'torder', 'Family' = 'tfamily',
            'Genus' = 'tgenus', 'Species' = 'tspecies', 'Strain' = 'tstrain'
          )
        )
      ),
      column(
        6,
        selectInput(
          'protein_rank', 'Protein rank', 
          c(
            'Superfamily' = 'psuperfamily', 'Family' = 'pfamily', 
            'Class' = 'pclass', 'Subclass' = 'psubclass', 'Group' = 'pgroup'
          )
        )
      )
    ),
      
    # Show a plot of the generated distribution
    fluidRow(
        tableOutput("matrix")
    ),
    fluidRow(
        plotOutput("summaryplot")
    )
  )
)
