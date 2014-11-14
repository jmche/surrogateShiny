
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggvis)
source("helper.R")

shinyUI(fluidPage(
  titlePanel("Surrogate Mutation Explorer for Breast Cancer Cell Lines"),

  tags$p(helpPopup(title="Surrogate Explorer Help",content = "<p>The control box is draggable if it gets in the way.
                          <p>Use the slider to adjust the alpha cutoff. <p>
                          Mouse over a box in the p-value matrix to show a tooltip with information of the surrogate
                          mutation. <p>
                          To see the surrogate mutation represented as a network, click on the box. <p>
                          To see if a gene has an alteration (mutation or copy number), click the checkbox.<p>
                          To sort by columns or rows, use the appropriate select box.  The clustering order for both rows and columns 
                          was done on pvalues filtered at an alpha of 0.05. Columns can also be sorted by number of alterations within
                          the sample.
                          (both mutations and copy number alterations). Rows also can be sorted by degree, or number of connections
                          a node has from the HPRD network.<p>
                   For more info, mail Ted Laderas (laderast@ohsu.edu).", trigger="click", placement="bottom"),
         "Code is available at",  tags$a(href="https://github.com/laderast/surrogateShiny", "github.com/laderast/surrogateShiny")),
  
  #sidebarLayout(
    absolutePanel(id = "controls", 
                  class = "modal", 
                  fixed = TRUE, draggable = TRUE, 
                  #position= "center",
                  top = 280, left="auto", right = 40, bottom = "auto", style="opacity: 0.8",
                  width = 300, height = "auto",
                div(
                h4("Controls"),
                
                sliderInput("pvalSlider", label =h5("Alpha"),
                  min = 0.0003, max = 0.1, value = 0.05),
                checkboxInput("mutBox", label="Show Genes With Mutations", value=FALSE),
                selectInput("selectCols", label = h5("Order columns"), 
                            choices = list("Alphabetical" = 1, "Clustering (pval < 0.05)" = 2,
                                           "Total Alterations (decreasing)" = 3), selected = 1),
                selectInput("selectRows", label = h5("Order rows"), 
                            choices = list("Alphabetical" = 1, "Clustering (pval < 0.05)" = 2,
                                           "Node Degree" = 3), 
                            selected = 1)
                
                #checkboxInput("clusterColumns", label=h4("Cluster Columns"), value=FALSE),
                
                
                )),
      absolutePanel(top = 100, left = 0, right = 0, fixed = FALSE,
        ggvisOutput("pvalPlot")
      )
    
  )
)

