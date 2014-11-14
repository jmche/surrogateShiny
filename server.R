
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggvis)
library(reshape2)
surrogateTable <- read.delim("surrogate-table.txt")
noGenes <- length(table(surrogateTable$Gene))
noSamples <- length(table(surrogateTable$Sample))

#load in precomputed column orders
levRow <- read.delim("rowOrder.txt", stringsAsFactors=FALSE)
levCol <- read.delim("colOrder.txt", stringsAsFactors=FALSE)

shinyServer(function(input, output, session) {

  ##surrTable is a reactive version of surrogateTable
  surrTable <- reactive({
    #input$pvalSlider
    pvalFilt <- ifelse(surrogateTable$pvalue < input$pvalSlider, 1, 0)
    if(input$mutBox){
      #add isMutated score x 2 to matrix
      #0 = nothing, 1 = surrogate mutation, 2 = mutated, 
      #3 = mutated + surrogate
      pvalFilt <- pvalFilt + (surrogateTable$isMutated *2)
      }
    pvalFilt <- factor(pvalFilt)
    return(data.frame(ID=surrogateTable$ID, Sample=surrogateTable$Sample, Gene=surrogateTable$Gene, pvalFilt = pvalFilt))
  })
  
  surrGene_tooltip <- function(x){
    if(is.null(x)) return(NULL)
    if(is.null(x$ID)) return(NULL)
    
    g <- surrogateTable[surrogateTable$ID == x$ID,]
    
    paste0("<b>Gene: </b>", g$Gene, "<br>",
           "<b>Sample: </b>",g$Sample, "<br>",
           "<b>P-value: </b>", g$pvalue, "<br>",
           "<b>Mutated Neighbors: </b>", g$neighbor, "<br>",
           "<b>All Neighbors: </b>", g$degree)
  }
  
  surrGraph_tooltip <- function(x){
    if(is.null(x)) return(NULL)
    if(is.null(x$ID)) return(NULL)
    
    img <- paste(x$ID,".svg", sep="")
    #select the appropriate row of surrogateTable
    g <- surrogateTable[surrogateTable$ID == x$ID,]
    
    #print(g)
    if(g[1,]$neighbor > 0){ 
      if(g[1,]$neighbor==1 & g[1,]$isMutated ==1){
        sprintf("%s is mutated in %s, but has no mutated neighbors.", g$Gene, g$Sample)
      }
      else{sprintf("<img src='%s' alt = '%s' width = 500 height = 500></img><br>
                   <img src='legend1.jpg' alt='legend' width=500></img>", img, x$ID)
      }  
       
    }
    else{sprintf("%s has no neighboring mutations for %s.", g$Gene, g$Sample)}
  }
  
  surrPlot <- reactive({
    #select column and row ordering from select boxes
    domX <- levCol[,as.numeric(input$selectCols)]
    domY <- levRow[,as.numeric(input$selectRows)]
    
    #rs <- tapply(surrTable()$pvalFilt, surrTable()$Gene, sum)

    #if(input$filterBox){
    #domY <- domY[domY %in% names(rs[rs>0])]
    #}
    if(input$mutBox){
      #specify the four categories if checkbox
      domFill <- c(0,1,2,3)
      rangeFill <- c("#132B43", "#56B1F7", 
                           "#FFC0CB", "#663399")
      scaleLabels <- c("Not Mutated", "Surrogate (p < alpha)", "Node Mutated", 
                       "Node Mutated + Surrogate (p < alpha)")
    }

    #otherwise, just show two categories
    else{domFill <- c(0,1)
      rangeFill <- c("#132B43", "#56B1F7")
      scaleLabels <- c("Not Mutated", "Surrogate (p < alpha)")
    }
    
    surrTable %>% ggvis(~Sample, ~Gene, fill=~pvalFilt) %>%
      #make the pvalue heatmap by calling layer_rects
      #note that each box has a key that corresponds to ID
      #in surrTable
      layer_rects(width = band(), height = band(), key := ~ID) %>%
      #display information when hovering
      add_tooltip(surrGene_tooltip, "hover") %>%
      #when clicking, display precomputed network graph
      add_tooltip(surrGraph_tooltip, "click") %>%
      #define the sample axis
      add_axis("x", properties = axis_props(labels = list(angle = 90)), orient="top", 
               title_offset = 80, tick_padding=30) %>%
      #define the gene axis
      add_axis("y", orient="left", title_offset = 80) %>%
      #define the fill scale
      scale_ordinal("fill", domain=domFill, range=rangeFill) %>%
      #define which samples to show and in what order using domX
      scale_ordinal("x", padding = 0, points = FALSE, domain = domX) %>%
      #define which genes to show and in what order using domY
      scale_ordinal("y", padding = 0, points = FALSE, domain = domY) %>%
      #add and label 
      add_legend("fill", title="Surrogate values", values=scaleLabels) %>%
      set_options(width= 5 * noGenes, height= 10 * noGenes)
    
  })
  
  
  #bind to shiny plot
  surrPlot %>% bind_shiny("pvalPlot") 
  
  
})
