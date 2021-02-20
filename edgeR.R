source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

library(edgeR)

exprsFile<-"/home/hema/Desktop/matrix2.csv"
#class<-as.vector(read.table(exprsFile), mode = "any")
exprs <- as.matrix(read.csv(exprsFile, header=TRUE, sep = ",",
                              row.names = 1,
                              as.is=TRUE))

exprsFile1<-"/home/hema/Desktop/m1.csv"
#class<-as.vector(read.csv(exprsFile1), mode = "any")
classcond <- as.data.frame(read.csv(exprsFile1, header=TRUE, sep = ",",
                              #row.names = 1,
                              as.is=TRUE))


edgeR.dgelist=DGEList(counts=exprs,  group = classcond$Condition)
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH") 



tstout = cbind(rownames(exprs),edgeR.adjpvalues)
write.table(tstout, "/home/hema/Desktop/matrix_output5.txt", sep = ",", col.names = F, row.names = F)

