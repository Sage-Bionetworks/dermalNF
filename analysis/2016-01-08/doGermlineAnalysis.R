source("../../bin/WGSData.R")


source("../../bin/dermalNFData.R")

##now re-annotate based on synapse id.
##see if we can find missing WGS
mapping=synTableQuery('SELECT "Patient","DnaID","WGS" FROM syn5556216')@values
mapping=mapping[which(!is.na(mapping$WGS)),]


samps<-synapseQuery("select * from entity where parentId=='syn5522788'")

annotes<-sapply(samps$entity.id,function(x) synGet(x,downloadFile=F)@annotations)
names(annotes)<-samp$entity.id

##now i did some updates of the sequence analysis
res1=getMutationStatsForGene("NF1")

res2=getMutationStatsForGene("TP53")
res3=getMutationStatsForGene("MUC6")
