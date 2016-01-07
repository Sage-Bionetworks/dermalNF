source("../../bin/WGSData.R")
#executed this command already
#storeMutationFiles()

library(data.table)

##now what?
allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
#try to create matrix.
##first read in all files
allmuts<-lapply(allsoms$entity.id,function(x) read.table(synGet(x)@filePath,sep=' ',header=T,quote='"'))

##now get genes
allgenes<-c()
for(a in allmuts)
  allgenes=union(allgenes,a$Hugo_Symbol)

#which gene is a mutant in which sample?
membership=sapply(allgenes,function(x) sapply(allmuts,function(y) x%in%y$Hugo_Symbol))

#which gene is a mutant in at least 2 sampls
res=apply(membership,2,function(x) length(which(x))==nrow(membership))
  
##get the top offenders? 

getMutationStatsForGene("NF1")
getMutationStatsForGene("MUC6")
getMutationStatsForGene("SMAD5")
getMutationStatsForGene("MAML3")