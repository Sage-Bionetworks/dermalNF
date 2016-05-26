##goal is to grab proteomics ratios and run them using the forest script ##from omics integrator

source("../../bin/dermalNFData.R")

#get proteomics data and annotations
prot.annotes<-protein_annotations()
prots <- prot_normalized()

##summarize counts in single forest, and in patient-specific forest.

prot.vals<-sapply(prots[,1],function(x) colSums(prots[which(prots[,1]==x),-1],na.rm=T))

prot.sums<-colSums(prot.vals)
prot.means<-colMeans(prot.vals)
prot.rank<-colSums(prot.vals)/sum(prot.vals)

parid='syn5804586'

write.table(data.frame(names(prot.sums),prot.sums),row.names=F,col.names=F,file='sumOfProteinsAcrossAllSamples.tab',sep='\t',quote=F)
write.table(data.frame(names(prot.means),prot.means),row.names=F,col.names=F,file='meanOfProteinsAcrossAllSamples.tab',sep='\t',quote=F)
write.table(data.frame(names(prot.rank),prot.rank),row.names=F,col.names=F,file='fracOfProteinsAcrossAllSamples.tab',sep='\t',quote=F)

allfs=c('sumOfProteinsAcrossAllSamples.tab','meanOfProteinsAcrossAllSamples.tab','fracOfProteinsAcrossAllSamples.tab')
for(f in allfs){
  synStore(File(f,parentId=parid),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-18/formatForestData.R',wasExecuted=TRUE)))
}
##now can we select fraction of samples by patient? 

sapply(unique(prot.annotes$patientId),function(x){
  if(x=='NULL')
    return(NULL)
  pid=gsub('CT0+','',x)
  synids=prot.annotes$synapseId[grep(paste('00',pid,'$',sep=''),prot.annotes$patientId)]
  if(length(synids)<2)
    return(NULL)
  pavals=prot.vals[synids,]
  pavals=pavals[,which(colSums(pavals,na.rm=T)>0)]
  prot.sums<-colSums(pavals)
  prot.means<-colMeans(pavals)
  prot.rank<-colSums(pavals)/sum(pavals)
  
  write.table(data.frame(names(prot.sums),prot.sums),row.names=F,col.names=F,file=paste('sumOfProteinsForPatient',pid,'.tab',sep=''),sep='\t',quote=F)
  write.table(data.frame(names(prot.means),prot.means),row.names=F,col.names=F,file=paste('meanOfProteinsForPatient',pid,'.tab',sep=''),sep='\t',quote=F)
  write.table(data.frame(names(prot.rank),prot.rank),row.names=F,col.names=F,file=paste('fracOfProteinsForPatient',pid,'.tab',sep=''),sep='\t',quote=F)
  
  allfs=c(paste('sumOfProteinsForPatient',pid,'.tab',sep=''),paste('meanOfProteinsForPatient',pid,'.tab',sep=''),paste('fracOfProteinsForPatient',pid,'.tab',sep=''))
  el=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-18/formatForestData.R'))
  ul=lapply(unique(synids),function(x) list(entity=x,wasExecuted=F))
    for(f in allfs){
      
    synStore(File(f,parentId=parid),
             used=ul,executed=el)
  }
  
})

