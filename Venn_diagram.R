install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.quad.venn(area1 = 3872, area2 = 11399, area3 = 10607, area4= 10238, n12 = 3862, n13 = 3630, n14 = 3666,
                n23 = 8736, n24 = 8449, n34 = 9855, n123 = 3622, n124 = 3657, n134 = 3595, n234 = 8274, 
                n1234 = 3587, category = c("DESeq", "edgeR", "Voom", "Vst"), lty = "blank", 
                fill = c("skyblue", "pink", "yellow","green"))

draw.quad.venn(area1 =315, area2 = 495, area3 = 330, area4 = 370, n12 = 145, n13 = 85, n14 = 105,
               n23 = 130, n24 = 95, n34 = 75, n123 = 50, n124 = 45, n134 = 30, n234 = 35, 
               n1234 = 5, category = c("DESeq", "edgeR", "Voom", "Vst"), lty = "blank", 
               fill = c("skyblue", "pink1", "yellow","red"))

## modified slightly from the example given in the documentation
## Example using a list of item names belonging to the
## specified group.
##
require(gplots) 
## construct some fake gene names..
oneName <- function() paste(sample(LETTERS,4,replace=TRUE),collapse="")
geneNames <- replicate(20502, oneName())

## 
GroupA <- sample(geneNames, 3872, replace=FALSE)
GroupB <- sample(geneNames, 11399, replace=FALSE)
GroupC <- sample(geneNames, 10610, replace=FALSE)
GroupD <- sample(geneNames, 10238, replace=FALSE)

venn(list(DESeq=GroupA,edgeR=GroupB,Voom=GroupC,Vst=GroupD))

plotAnimals(c("d", "c", "l", "s"), category = c("Dog People", "Cat People", 
                                                "Lizard People", "Snake People"), lty = "blank", fill = c("skyblue", "pink1", 
                                                                                                          "mediumorchid", "orange"))
