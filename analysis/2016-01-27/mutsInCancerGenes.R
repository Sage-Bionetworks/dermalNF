##compare cancer genes for mutational data.


source("../../bin/WGSData.R")

##now store the file
#cfile='../../data/Census_allTue Jan 19 18-58-56 2016.csv'
#synStore(File(cfile,parentId='syn4984723'))


##now get
cancer.genes=read.table(synGet('syn5610781')@filePath,sep=',')
allstats<-sapply(cancer.genes[,1],getMutationStatsForGene,doPlot=FALSE)