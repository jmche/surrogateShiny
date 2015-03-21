library(shiny)
library(ggvis)
library(Rgraphviz)
library(graph)
library(gplots)
if(file.exists("www/netgraph.svg")){file.remove("www/netgraph.svg")}
source("helper.R")
source("plot-results.R")

shinyUI(fluidPage(
  titlePanel("Surrogate Mutation Explorer for Breast Cancer Cell Lines"),

  tags$p(
    "Code is available at",  tags$a(href="https://github.com/laderast/surrogateShiny", "github.com/laderast/surrogateShiny", 
                                    tags$br(),
          helpPopup(title="Surrogate Explorer Help",content = "<p>The control box is draggable if it gets in the way.
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
                   For more info, mail Ted Laderas (laderast@ohsu.edu).", trigger="click", placement="bottom") 
         )),
  absolutePanel(top = 100, left = 0, right = 0, fixed = FALSE,
                ggvisOutput("pvalPlot")
  ),

  absolutePanel(id = "netgraph", 
                class = "panel panel-body", 
                fixed = FALSE, 
                #draggable = TRUE, 
                #position= "center",
                top = 280, left=765, right ="auto", bottom = "auto", 
                #style="opacity: 0.8",
                width = 500, height = "auto",
                
                #if(file.exists("www/netgraph.svg")){
                h5("Surrogate Graph (Click to Enlarge)"),
                a(imageOutput("imageNet"), href="graphView.html", 
                  target="_blank"),  
                br(),
                img(src='legend1.jpg', width=400)
                #}
  ),
  
  #sidebarLayout(
    absolutePanel(id = "controls", 
                  class = "panel panel-body", 
                  fixed = TRUE, 
                  draggable = TRUE, 
                  #position= "center",
                  top = 165, left="auto", right = 10, bottom = "auto", 
                  style="opacity: 0.8; background-color: lightgrey",
                  width = 250, height = "auto",
                  #tags$style("body {background-color: grey;}"),
                  
                  h5("Controls (Draggable)"),
                  
                  sliderInput("pvalSlider", label =h6("Alpha"),
                              min = 0.0001, max = 0.1, value = 0.05, step = 0.01),
                  checkboxInput("mutBox", label="Show Genes With Mutations", value=FALSE),
                  selectInput("selectCols", label = h6("Order columns"), 
                              choices = list("Alphabetical" = 1, "Clustering (pval < 0.05)" = 2,
                                             "Total Alterations (decreasing)" = 3), selected = 2),
                  selectInput("selectRows", label = h6("Order rows"), 
                              choices = list("Alphabetical" = 1, "Clustering (pval < 0.05)" = 2,
                                             "Node Degree" = 3), 
                              selected = 2)
                  
                
                #checkboxInput("clusterColumns", label=h4("Cluster Columns"), value=FALSE),
                
                
                )
    
  )
)

