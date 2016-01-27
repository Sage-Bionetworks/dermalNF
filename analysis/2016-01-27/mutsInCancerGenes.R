##compare cancer genes for mutational data.


source("../../bin/WGSData.R")

##now store the file
#cfile='../../data/Census_allTue Jan 19 18-58-56 2016.csv'
#synStore(File(cfile,parentId='syn4984723'))

require(parallel)
##now get
cancer.genes=read.table(synGet('syn5610781')@filePath,sep=',')

impact=c("LOW",'MODERATE','HIGH')
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  print(paste('Selecting from',nrow(allsoms),'mutation files'))
  allsoms=allsoms[unlist(sapply(impact,grep,allsoms$entity.name)),]
  print(paste("Found",nrow(allsoms),'with',paste(impact,collapse=' or '),'impact'))
  som.germ=getAllMutData(allsoms)

allstats<-mclapply(as.character(cancer.genes[,1]),function(x) try(getMutationStatsForGene(gene=g,doPlot=FALSE,som.germ=som.germ)),mc.cores=24)
