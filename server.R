library(shiny)
library(ggvis)
library(dplyr)
#library(reshape2)
library(Rgraphviz)
load("tcgaResult.RData")
load("intome.RData")
#cutoffMat <- read.delim("cutoff-matrix.txt")
surrogateTable <- read.delim("surrogate-table.txt")
subTypeMelt <- read.delim("subtype-melted.txt")
#subTypeMelt$pvalue <- subTypeMelt$pvalue * -3
surrogateTable <- rbind(surrogateTable, subTypeMelt)
geneIntTable <- read.delim("geneIntTable.txt")
noGenes <- length(table(surrogateTable$Gene))
noSamples <- length(table(surrogateTable$Sample))

#load in precomputed column orders
levRow <- read.delim("rowOrder.txt", stringsAsFactors=FALSE)
levCol <- read.delim("colOrder.txt", stringsAsFactors=FALSE)

shinyServer(function(input, output, session) {

  plotNetGraphsFromID <- function(ID, gisticCopyCalls=NULL, surrogateResults, 
                                  intome, fileOut) {
    #source("~/Code/surrogate/plot-results.R")
    NodeName <- as.character(surrogateTable[surrogateTable$ID==ID,"NodeName"])
    Sample <- as.character(surrogateTable[surrogateTable$ID==ID,"Sample"])
    GeneName <- as.character(surrogateTable[surrogateTable$ID==ID,"Gene"])
    #print(paste(NodeName,Sample))
    renderSubNetSimple(NodeName, Sample, GeneName, intome, gisticCopyCalls, 
                       surrogateResults, geneIntTable, fileOut)
  }

  surrGene_tooltip <- function(x){
    if(is.null(x)) return(NULL)
    if(is.null(x$ID)) return(NULL)
    if(x$Gene == "Subtype") return(NULL)
    
    g <- surrogateTable[surrogateTable$ID == x$ID,]
    
    paste0("<b>Click to Graph in Sidepanel</b><br>",
           "<b>Gene: </b>", g$Gene, "<br>",
           "<b>Sample: </b>",g$Sample, "<br>",
           "<b>P-value: </b>", g$pvalue, "<br>",
           "<b>Mutated Neighbors: </b>", g$neighbor, "<br>",
           "<b>All Neighbors: </b>", g$degree        
           )
  }

  
  surrGraphDynamic_tooltip <- function(x){
    if(is.null(x)) return(NULL)
    if(is.null(x$ID)) return(NULL)
    if(x$Gene == "Subtype") return(NULL)
    
    #img <- paste(x$ID,".svg", sep="")
    #select the appropriate row of surrogateTable
    g <- surrogateTable[surrogateTable$ID == as.character(x$ID),]
    
    g <- g[1,]
    if(g[1,]$neighbor > 0){ 
      if(g[1,]$neighbor==1 & g[1,]$isMutated ==1){
        sprintf("%s is mutated in %s, but has no mutated neighbors.", g$Gene, g$Sample)
      }
      else{
         
          outfile <- "netgraph.svg"
          fileOut <- paste("www/", outfile, sep="")
          g <- surrogateTable[surrogateTable$ID == as.character(x$ID),]
          #print(g$ID)
          try(plotNetGraphsFromID(as.character(g$ID), surrogateResults=tcgaResult, intome=intome,fileOut=fileOut), silent=TRUE)
          #plotNetGraphsFromID(as.character(g$ID), surrogateResults=tcgaResult, intome=intome,fileOut=fileOut)
          
          paste0("<b>Graph has been generated in sidepanel.</b><br>",
            "<b>Gene: </b>", g$Gene, "<br>",
                 "<b>Sample: </b>",g$Sample, "<br>",
                 "<b>P-value: </b>", g$pvalue, "<br>",
                 "<b>Mutated Neighbors: </b>", g$neighbor, "<br>",
                 "<b>All Neighbors: </b>", g$degree
                 )
          
     
          
          #sprintf("<img src='%s' alt = '%s' width = 500 height = 500></img><br>
          #         <img src='legend1.jpg' alt='legend' width=500></img>", outfile, g$ID)
      }  
      
    }
    else{sprintf("%s has no neighboring mutations for %s.", g$Gene, g$Sample)}
  }
  
  
  
surrPlot <- reactive({
    #select column and row ordering from select boxes
    domX <- levCol[,as.numeric(input$selectCols)]
    domY <- levRow[,as.numeric(input$selectRows)]
    
    subTypeVals <- c(-1,-2,-3,-10) 
    subTypeCols <- c("FF3333", "66CC00", "FFCC00", "000000")
    subTypeLabs <- c("Basal A", "Basal B", "Luminal", "No Subtype")
    
    #rs <- tapply(surrTable()$pvalFilt, surrTable()$Gene, sum)
    domFill <- c(0,1)
    rangeFill <- c("#132B43", "#56B1F7")
    scaleLabels <- c("Not Mutated", "Surrogate (p < alpha)")
    
    if(input$mutBox){
      domFill <- c(domFill, c(2, 3))
      rangeFill <-  c(rangeFill, c("#FFC0CB", "#663399"))
      scaleLabels <- c(scaleLabels, c("Node Mutated", "Node Mutated + Surrogate (p < alpha)"))
    }

    domFill <- c(domFill, subTypeVals)
    rangeFill <- c(rangeFill, subTypeCols)
    scaleLabels <- c(scaleLabels, subTypeLabs)


#})
    
    
    pvalFiltFun <-  function(pvalue, isMutated) {
      ifelse(pvalue < input$pvalSlider, 1,0) + (input$mutBox * isMutated*2)
    }  

    surrogateTable %>% 
      select(ID, Gene, Sample, pvalue, isMutated) %>% 
      #pvalFilt is a factor ranging from 0-3
      #use negative values to denote subtype
      mutate(pvalFilt = factor(ifelse(pvalue < 0, pvalue, pvalFiltFun(pvalue, isMutated)))) %>%
      ggvis(~Sample, ~Gene, fill=~pvalFilt) %>%
      #make the pvalue heatmap by calling layer_rects
      #note that each box has a key that corresponds to ID
      #in surrTable
      layer_rects(width = band(), height = band(), key := ~ID) %>%
                  #fill.brush := "green") %>%
      #display information when hovering
      add_tooltip(surrGene_tooltip, "hover") %>%
      #when clicking, display precomputed network graph
      add_tooltip(surrGraphDynamic_tooltip, "click") %>%    
      #add_tooltip(surrGraph_tooltip, "click") %>%
      
      #brush functionality for summarizing networks
      #handle_brush(function(items, page_loc, session, ...){
      #   print(items$ID)
      #   val <- summarizeNetworks(items$ID, intome, tcgaResult$mutCopyFrames, surrogateTable, geneIntTable)
      #   print(val)
         #if(!is.null(val)){
         #show_tooltip(session, html="<img src='netgraph.svg' alt = 'summary graph' width = 500 height = 500></img>")}
      #}) %>%

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
      add_legend("fill", title="Legend", values=scaleLabels) %>%
      #add_legend("fill", title="PAM50 Subtypes", values=subTypeLabs) %>%
      set_options(width= 15 * (noSamples), height= 10 * (noGenes +1)) 
    
 })
  
  #bind to shiny plot
  surrPlot %>% bind_shiny("pvalPlot") 

  #poll www/netgraph.svg to see if the netgraph is updated
  svgGraph <- reactiveFileReader(1000, session, 'www/netgraph.svg', 
                                 function(x) {return("www/netgraph.svg")})
  
  output$imageNet <- renderImage({
    list(src = svgGraph(),
         contentType = "image/svg+xml", width= 500, height=500
         ) 
  },deleteFile=FALSE)
  
  
  
})
