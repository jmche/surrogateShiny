source("plot-results.R")
load("tcgaResult.RData")
surrogateTable <- read.delim("surrogate-table.txt")

testVec <- c("BT20-BRCA1" , "BT20-BRF1"  , "BT20-CASP2" , "BT20-CASP3",  
             "BT474-BRCA1", "BT474-BRF1",  "BT474-CASP2", "BT474-CASP3",
             "BT483-BRCA1", "BT483-BRF1",  "BT483-CASP2", "BT483-CASP3" ,"BT549-BRCA1", "BT549-BRF1" ,
 "BT549-CASP2", "BT549-CASP3")

testVec <- c("HCC1395-APLP2","HCC1395-APPL","HCC1395-AR","HCC1395-ATN1","HCC1395-ATP2B2","HCC1428-APLP2", 
"HCC1428-APPL",   "HCC1428-AR","HCC1428-ATN1",   "HCC1428-ATP2B2", "HCC1569-APLP2",  "HCC1569-APPL",  
"HCC1569-AR",     "HCC1569-ATN1",   "HCC1569-ATP2B2", "HCC1806-APLP2",  "HCC1806-APPL",   "HCC1806-AR",    
"HCC1806-ATN1",   "HCC1806-ATP2B2", "HCC1937-APLP2",  "HCC1937-APPL",   "HCC1937-AR",     "HCC1937-ATN1",  
"HCC1937-ATP2B2", "HCC1954-APLP2",  "HCC1954-APPL",   "HCC1954-AR",     "HCC1954-ATN1",   "HCC1954-ATP2B2",
"HCC202-APLP2",   "HCC202-APPL",    "HCC202-AR",      "HCC202-ATN1",    "HCC202-ATP2B2" )

summarizeNetworks(testVec, intome, tcgaResult$mutCopyFrames, 
                  surrogateTable, geneIntTable, NULL)

testVec <- as.character(surrogateTable[surrogateTable$Gene == "ESR1","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "BRCA1","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "RB1","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "MAPK9","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "SMARCB1","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "PTK2","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "GSK3B","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "AR","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "UBB","ID"])
testVec <- as.character(surrogateTable[surrogateTable$Gene == "CTNNB1","ID"])




testVec <- as.character(surrogateTable[surrogateTable$Sample == "SKBR3" & surrogateTable$pvalue < 0.05,"ID"])
testVec <- as.character(surrogateTable[surrogateTable$Sample == "HCC202" & surrogateTable$pvalue < 0.05,"ID"])
testVec <- as.character(surrogateTable[surrogateTable$Sample == "BT20" & surrogateTable$pvalue < 0.01,"ID"])




summarizeNetworks(testVec, intome, tcgaResult$mutCopyFrames, surrogateTable, geneIntTable, NULL)

