##get original proteomics filenames

pd='/scratch/DERMALNF/disk1/proteomics_data/'
all.files=list.files(pd)

require(synapseClient)
synapseLogin()
syn.files= synapseQuery('SELECT name,ID,patientID,tissueID FROM entity WHERE parentId=="syn4984949"')

#map sample names to files
for(f in all.files){
  tab<-read.table(paste(pd,f,sep='/'),sep = '\t',header = T)
  samps=sapply(colnames(tab)[6:8], function(x) gsub('.','_',x,fixed=T))
  idx=sapply(samps,grep,syn.files$entity.name)
  synids=syn.files$entity.id[idx]
  for(s in synids){
    sf=synGet(s,downloadFile=F)
    synSetAnnotations(sf) <- list(originalBatch=f,dataType='iTRAQ',tissueType='Tissue',
                                  patientID=syn.files$entity.patientID[match(s,syn.files$entity.id)],
                                  tissueID=syn.files$entity.tissueID[match(s,syn.files$entity.id)])
    synStore(sf)
  }
}
