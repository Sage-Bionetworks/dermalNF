##get original proteomics filenames

all.files=list.files('/scratch/DERMALNF/disk1/proteomics_data/')

require(synapseClient)

syn.files= synapseQuery('SELECT name,ID,patientID,tissueID FROM entity WHERE parentId=="syn4984949"')

#map sample names to files
for(f in all.files){
  tab<-read.table(f,sep = '\t',header = T)
  samps=colnames(tab)[6:8]
  
}
