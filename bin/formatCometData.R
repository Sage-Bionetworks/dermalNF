###goal is to download dermal somatic mutations and run networks

require(synapseClient)

synapseLogin()
require(dplyr)

##get all cancer mutations
#allCancerMuts=synGet('syn5611520')

allMuts <- synGet("syn5839666")
##get all mutations

##get expressed genes!
source("../../bin/dermalNFData.R")
rna.ex=rna_count_matrix(minCount=3)
expressed.genes=rownames(rna.ex)
print(paste('Found',length(expressed.genes),'and reducing mutations accordingly'))


getCountsFromTable<-function(synfile,mutationType=c('germline','somatic'),prefix='',includeSilent=FALSE){
  tab <- read.table(synfile@filePath,sep = '\t',header = T)
  
  if(!includeSilent)
    tab <- tab[which(tab$Mutation_Type!='Silent'),]
  
  mut.idx = which(sapply(tab$Mutation_Status,tolower)%in%mutationType)
  tab <- tab[mut.idx,]

  tab<- tab[which(tab$Hugo_Symbol%in%expressed.genes),]
  
  mut.genes = tab %>% group_by(Sample_ID) %>% summarize(Genes = distinct(Hugo_Symbol))
  
  fname= paste(prefix,paste(mutationType,collapse='_and_'),ifelse(includeSilent,'withSilent','nonSilent'),'mutGenesPerPatient.tab',sep='_')
  write.table(mut.genes,file=fname,col.names=F,row.names=F,quote=F)
}

#for(mutType in c('germline','somatic',c('germline','somatic'))){
mutType<-'germline'
  for(silent in c(TRUE,FALSE)){
#    getCountsFromTable(allCancerMuts,mutationType=mutType,prefix='allCancerGenes',includeSilent = silent)
    getCountsFromTable(allMuts,mutationType=mutType,prefix='allGenes',includeSilent = silent)
    
  }

