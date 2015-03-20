library(gplots)
library(graph)
library(Rgraphviz)

renderSubNet <- function(NodeName, intome=intome, cellLine, 
                         fileprefix="Subnetwork", gisticCopyCalls=gisticCopyCalls, 
                         mutCopyFrames=mutCopyFrames, RPPANodes=RPPANodes, 
                         cellResults=cellResults, mappedDoms=mappedDoms,
			 proVEANgenemerge){

  library(Rgraphviz)
  library(BioNet)
  
  cellFrame <- gisticCopyCalls[,c("GeneName", "NodeName", cellLine)]

  mutCopyFrame <- mutCopyFrames[[cellLine]]
  
  mutDoms <- mappedDoms[mappedDoms$CellLine==cellLine,]
  
  proveanRes <- proVEANgenemerge[proVEANgenemerge$CellLine==cellLine,]
  
  mutDomMerge <- merge(mutCopyFrame, mutDoms, by.x="Gene", by.y="HGNC_ID")
    
  mutDomMerge <- mutDomMerge[,c("Gene", "NodeName", "hmm_name", "clan")]
  
  mutDomMerge <- unique(mutDomMerge)
  
  #merge filtered results
  pvmuts <- merge(mutCopyFrame, proveanRes, by.x="Gene", by.y="HGNC_ID")
  pvmuts <- pvmuts[,c("NodeName","PREDICTION..cutoff.0.05.")]
  pvmuts <- unique(pvmuts)

  nodenet <- c(NodeName,intersect(as.character(inEdges(NodeName, intome)[[1]]), 
                                  as.character(mutCopyFrame$NodeName)))  
  
  #grab those nodes with mutations and those nodes with copy alterations
  brcamut <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$mutations!=0,"NodeName"]))
  #brcacopy <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$copy!=0,"NodeName"]))
  
  #pull labels
  labeltab <- geneIntTable[geneIntTable$NodeName %in% nodenet,]
  labels <- as.character(labeltab$Gene)
  names(labels) <- labeltab$NodeName
  
  for(nod in mutDomMerge$NodeName){
     annot <- mutDomMerge[mutDomMerge$NodeName == nod,]
     #print(annot)
     if(dim(annot)[1] > 0){
     labels[nod] <- paste(labels[nod], "\n", paste(annot$hmm_name), "\n" , paste(annot$clan), sep="")
     }
  }
  
  nodenetwork <- subNetwork(nodenet, intome)
  
  nodeRenderInfo(nodenetwork) <- list(shape="circle", 
	iheight=.5, iwidth= .5, 
	fixedsize=FALSE, label = as.list(labels))
  inRPPA <- nodenet[nodenet %in% RPPANodes]
  
  #make those RPPA proteins a diamond
  RPPAshape <- list()
  for(RP in inRPPA){RPPAshape[[RP]] <- "diamond"}
  nodeRenderInfo(nodenetwork) <- list(shape=RPPAshape)
    
  #UACCgistframe <- data.frame(rownames(copyGistFiltered),copyGistFiltered$UACC812)
  #colnames(UACCgistframe) <- c("GeneName", "copyStatus")
  #UACCcopynodestatus <- merge(UACCgistframe, geneIntTable, by.x="GeneName", by.y="Gene")
  
  nodecopyframe <- cellFrame[as.character(cellFrame$NodeName) %in% nodenet,]
  rownames(nodecopyframe) <- nodecopyframe$NodeName
  colnames(nodecopyframe) <- c("GeneName", "NodeName", "copyStatus")
  
  nodecolors <- list()
  
  for(nd in rownames(nodecopyframe)){
    ndcol <- "white"
    if(nodecopyframe[nd, "copyStatus"] == 1){ndcol <- "pink"}
    if(nodecopyframe[nd, "copyStatus"] == -1){ndcol <- "lightgreen"}
    nodecolors[[nd]] <- ndcol  
  }
  
  for(nd in brcamut) {nodecolors[[nd]] <- "lightblue"}

  pvnodes <- intersect(as.character(pvmuts$NodeName), as.character(brcamut))	

  #annotate damaging mutations as a darker blue
  for(nod in pvnodes){
     pvscores <- as.character(pvmuts[pvmuts$NodeName == nod,"PREDICTION..cutoff.0.05."])
     print(pvscores)
     if(length(grep("Damaging", pvscores) > 0)){
	nodecolors[[nod]] <- "deepskyblue3"
     }  
  }

  #distMat <- distMatrices[[cellLine]]
  #colnames(distMat) <- RPPANodes

  cellRes <- cellResults[[cellLine]]

  NodDegree <- cellRes[cellRes$NodeName ==NodeName,"degree"]
  NeighborMuts <- cellRes[cellRes$NodeName == NodeName,"neighborVec"]
  pval <- cellRes[cellRes$NodeName == NodeName, "pvalue"]

  nodenetwork <- layoutGraph(nodenetwork)
  nodeRenderInfo(nodenetwork) <- list(fill = nodecolors)
  pdf(height=7, width=7, file=paste(cellLine,"-",fileprefix,"-",NodeName,".pdf", sep=""))
  #plot.new()
  layout(matrix(c(1,2),2,1), heights=c(2,1))
  plotTitle <- paste(NodeName, " (",cellLine," ",", p=", pval, ", n=", NeighborMuts, ", d=", NodDegree, ")", sep="")
  renderGraph(nodenetwork)
  hist(distMat[,NodeName], main=plotTitle, xlab="Neighbor Mutations", ylab = "Frequency")
  dev.off()

}

#cellResults <- cellResults[-8]

renderSubNetSimple <- function(NodeName, cellLine, GeneName, intome=intome, 
                         gisticCopyCalls=NULL, resultObj, geneIntTable,
                         fileOut=NULL){
  
  library(Rgraphviz)
  library(graph)
  
  #grab gistic status
  if(!is.null(gisticCopyCalls)){
  cellFrame <- gisticCopyCalls[,c("GeneName", "NodeName", cellLine)]
  }
  
  mutCopyFrame <- resultObj$mutCopyFrames[[cellLine]]
  
  cellResults <- resultObj$cellResults
  names(cellResults) <- make.names(names(cellResults))
  
  nodenet <- c(NodeName,intersect(as.character(inEdges(NodeName, intome)[[1]]), 
                                  as.character(mutCopyFrame$NodeName)))  
  
  #grab those nodes with mutations and those nodes with copy alterations
  brcamut <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$mutations!=0,"NodeName"]))
  #brcacopy <- intersect(nodenet, as.character(mutCopyFrame[mutCopyFrame$copy!=0,"NodeName"]))
  
  #pull labels
  labeltab <- geneIntTable[geneIntTable$NodeName %in% nodenet,]
  labels <- as.character(labeltab$Gene)
  names(labels) <- labeltab$NodeName
    
  nodenetwork <- subGraph(nodenet, intome)
  
  nodeRenderInfo(nodenetwork) <- list(shape="circle", 
                                      iheight=.5, iwidth= .5, 
                                      fixedsize=FALSE, label = as.list(labels))
  #inRPPA <- nodenet[nodenet %in% RPPANodes]
  
  #make Node protein a diamond
  RPPAshape <- list()
  RPPAshape[[NodeName]] <- "diamond"  
  #for(RP in inRPPA){RPPAshape[[RP]] <- "diamond"}
  nodeRenderInfo(nodenetwork) <- list(shape=RPPAshape)
  
  if(!is.null(gisticCopyCalls)){
  nodecopyframe <- cellFrame[as.character(cellFrame$NodeName) %in% nodenet,]
  rownames(nodecopyframe) <- nodecopyframe$NodeName
  colnames(nodecopyframe) <- c("GeneName", "NodeName", "copyStatus")
  }
  else{
    nodecopyframe <- mutCopyFrame[,c("Gene", "NodeName", "copyStatus")]
  }
  
  nodecolors <- list()
  
  for(nd in rownames(nodecopyframe)){
    ndcol <- "white"
    if(nodecopyframe[nd, "copyStatus"] > 0){ndcol <- "pink"}
    if(nodecopyframe[nd, "copyStatus"] < 0){ndcol <- "lightgreen"}
    nodecolors[[nd]] <- ndcol  
  }
  
  for(nd in brcamut) {nodecolors[[nd]] <- "lightblue"}  
  
  cellRes <- cellResults[[cellLine]]
  
  NodDegree <- cellRes[cellRes$NodeName ==NodeName,"degree"]
  NeighborMuts <- cellRes[cellRes$NodeName == NodeName,"neighborVec"]
  pval <- cellRes[cellRes$NodeName == NodeName, "pvalue"]
  
  nodenetwork <- layoutGraph(nodenetwork)
  if(NeighborMuts > 0 & length(nodes(nodenetwork)) > 1){
    nodeRenderInfo(nodenetwork) <- list(fill = nodecolors)
  if(is.null(fileOut)){
    fileOut <- paste(cellLine,"-",GeneName,".svg", sep="")
  }
  svg(height=7, width=7, filename=fileOut)
  #plot.new()
  #layout(matrix(c(1,2),2,1), heights=c(2,1))
  renderGraph(nodenetwork)
  plotTitle <- paste(NodeName, " (",cellLine," ",", p=", pval, ", n=", NeighborMuts, ", d=", NodDegree, ")", sep="")
  title(plotTitle)
  #hist(distMat[,NodeName], main=plotTitle, xlab="Neighbor Mutations", ylab = "Frequency")
  dev.off()
  }  
}


summarizeNetworks <- function(IDvec, intome, mutCopyFrames, surrogateTable, geneIntTable, fileOut="www/netgraph.svg"){
  #library(gplots)
  if(length(IDvec) < 1){return(NULL)}
  
  nodeTable <- surrogateTable[surrogateTable$ID %in% IDvec,]
  
  fullNet <- NULL
  surrNodes <- unique(as.character(nodeTable$NodeName))
  #for each row of the table, grab neighbors of the surrogate node
  for(i in 1:nrow(nodeTable)){
    Sample <- nodeTable[i,"Sample"]
    NodeName <- as.character(nodeTable[i,"NodeName"])
    nodenet <- c(NodeName,intersect(as.character(inEdges(NodeName, intome)[[1]]), 
                                  as.character(mutCopyFrames[[Sample]]$NodeName)))
    fullNet <- c(fullNet, nodenet)
    }
  
  #count the number of times a node is mutated
  fullNet <- table(fullNet)
  
  surrNodeCount <- table(as.character(nodeTable$NodeName)) +1
  
  #subtract 1 from the count for the surrogate nodes, since we already count them
  fullNet[names(surrNodeCount)] <- fullNet[names(surrNodeCount)] - surrNodeCount 

  #don't graph if the network only has 1 node or less  
  if(sum(fullNet) < 2 | length(fullNet) < 2){return(NULL)}
    
  numColors <- max(fullNet) + 1
  #colorPanel <- topo.colors(numColors)
  colorPanel <- colorpanel(n=numColors, low="white", high="steelblue")
  
  #assign colors according to degres
  nodeColorList <- list()
  for(nd in names(fullNet)){
    nodeColorList[[nd]] <- colorPanel[fullNet[nd]]
  }
  
  #extract the relevant subgraph from the interactome
  fullGraph <- subGraph(names(fullNet), intome)
  
  labeltab <- geneIntTable[geneIntTable$NodeName %in% names(fullNet),]
  labels <- as.character(labeltab$Gene)
  names(labels) <- labeltab$NodeName
  
  #print(labels)
  
  #default node render properties
  nodeRenderInfo(fullGraph) <- list(shape="circle", 
                                      iheight=.5, iwidth= .5, 
                                      fixedsize=FALSE, label = as.list(labels))
  #inRPPA <- nodenet[nodenet %in% RPPANodes]
  
  #make surrogate node proteins a diamond
  shapeList <- list()
  for(surr in surrNodes){
    shapeList[[surr]] <- "diamond"  
  }
  #for(RP in inRPPA){RPPAshape[[RP]] <- "diamond"}
  nodeRenderInfo(fullGraph) <- list(shape=shapeList)
  
  #color each node by degree of mutation
  fullGraph <- layoutGraph(fullGraph, layoutType="fdp")
  
  nodeRenderInfo(fullGraph) <- list(fill = nodeColorList)
  if(!is.null(fileOut)){
    svg(height=7, width=7, filename = fileOut)
    renderGraph(fullGraph)
    title("Summary Network")
    dev.off()
  }
  else{renderGraph(fullGraph)
  }
}
