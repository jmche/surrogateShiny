load("tcgaResult.RData")
source("~/Code/surrogate/processSurrogateResults.R")
surrogateTable <- buildSurrogateTable(tcgaResult)
write.table(surrogateTable, "surrogate-table.txt", quote=F, sep="\t")

surrogateTable <- read.delim("surrogate-table.txt")

pvalMat <- dcast(surrogateTable, Gene ~ Sample, value.var="pvalue")
rownames(pvalMat) <- pvalMat$Gene
pvalMat <- pvalMat[,-1]
pvalFilt <- filterPvalMat(pvalMat)

clust <- heatmap(pvalFilt, scale="none")
rowClust <- clust$rowInd
colClust <- clust$colInd

mutCopyFrames <- tcgaResult$mutCopyFrames
numMuts <- unlist(lapply(mutCopyFrames, function(x){sum(x$combined)}))

#build column and row orders
#levRow: 1: alphabetic order (default) 2: clustered genes 3: Total degree of gene in HPRD
levRow <- data.frame(default=as.character(rownames(pvalFilt)),
                     hclust=as.character(rownames(pvalFilt)[rowClust]),
                     degree=names(sort(tapply(surrogateTable$degree, surrogateTable$Gene, mean), 
                     decreasing=TRUE)))
#levCol 1: alphabetic order (default) 2: clustered samples 3: Total # of mutations in sample
levCol <- data.frame(default=as.character(colnames(pvalFilt)), 
                     hclust=as.character(colnames(pvalFilt)[colClust]),
                     numberMutations=make.names(names(numMuts[order(numMuts,decreasing = TRUE)])))

write.table(levRow, "rowOrder.txt", sep="\t", quote=F)
write.table(levCol, "colOrder.txt", sep="\t", quote=F)
