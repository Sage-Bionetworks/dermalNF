source("../../bin/WGSData.R")

##now try to upload germline mutations
res=storeGermlineMutationFiles()

##now get germline info? 

source("../../bin/dermalNFData.R")
cnv_seg=cnv_segmented()

##now re-annotate based on synapse id.
mapping=synTableQuery('SELECT "Patient","DnaID","SNPArray" FROM syn5556216')@values
mapping=mapping[which(!is.na(mapping$SNPArray)),]
samps<-sapply(as.character(cnv_seg$ID),function(x) {
  idx=match(x,mapping$SNPArray)
  #print(idx)
  return(paste('Patient',mapping$Patient[idx],'Samp',gsub(' ','',mapping$DnaID[idx]),sep='_'))})

new_seg=cnv_seg
new_seg$ID=samps
write.table(new_seg,'CNV_Segmented_sd2.seg',row.names=F,quote=F,sep='\t')

##now i did some updates of the sequence analysis
res1=getMutationStatsForGene("NF1")

res2=getMutationStatsForGene("TP53")
res3=getMutationStatsForGene("MUC6")
