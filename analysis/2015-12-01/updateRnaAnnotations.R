##first goal: update feature counts files with proper annotations

library(synapseClient)
synapseLogin()

countfiles=synapseQuery("select id,name from entity where parentId=='syn4984701'")
countfiles=countfiles[grep('_featureCounts.txt',countfiles$entity.name),]

updated.files<-sapply(countfiles$entity.id,function(x){
  fi<-synGet(x,downloadFile = FALSE)
  #get provenance, bam file used, for each file
  
  allused=sapply(used(fi),function(y) synGet(y$reference$targetId,downloadFile=FALSE))
  fnames=sapply(allused,function(y) y@properties$name)
  #now which used file was actually the 
  bf=allused[[grep('.bam',fnames)]]
  #get annotations for bam file
  
  #now just re-assign the annotations
  annotations(fi)<-annotations(bf)
  ##re-upload annotations, keeping provenance!
  synStore(fi,used=used(fi))
  return(fi)
})

##NOW we can re-do the mapping table, so it will reference the count data
source('../../bin/crossDataMapping.R')
storePatientMappingDataInTable()


