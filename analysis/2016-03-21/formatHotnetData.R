###goal is to download dermal somatic mutations and run networks

require(synapseClient)

synapseLogin()
require(dplyr)

##get all cancer mutations
allCancerMuts=synGet('syn5611520')

allMuts <- synGet("syn5713423")
##get all mutations


getCountsFromTable<-function(synfile,mutationType=c('germline','somatic'),prefix='',includeSilent=FALSE){
  tab <- read.table(synfile@filePath,sep = '\t',header = T)
  
  if(!includeSilent)
    tab <- tab[which(tab$Mutation_Type!='Silent'),]
  
  mut.idx = which(sapply(tab$Mutation_Status,tolower)%in%mutationType)
  tab <- tab[mut.idx,]
  
  num.patients = tab %>% group_by(Hugo_Symbol) %>% summarize(NumPatients = n_distinct(Sample_ID))
  
  fname= paste(prefix,paste(mutationType,collapse='_and_'),ifelse(includeSilent,'withSilent','nonSilent'),'numMutsPerGene.tab',sep='_')
  write.table(num.patients,file=fname,col.names=F,row.names=F,quote=F)
}

for(mutType in c('germline','somatic',c('germline','somatic'))){
  for(silent in c(TRUE,FALSE)){
    getCountsFromTable(allCancerMuts,mutationType=mutType,prefix='allCancerGenes',includeSilent = silent)
    getCountsFromTable(allMuts,mutationType=mutType,prefix='allGenes',includeSilent = silent)
    
  }
}
