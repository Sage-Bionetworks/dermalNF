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
all.muts<-getAllMutData(impact='HIGH')

#res=storeSomMutationFiles(impact='MODERATE')
res=getMutationStatsForGene("OR4C5")
res=getMutationStatsForGene("FOXD1")
res=getMutationStatsForGene("MUC6")

res=getMutationStatsForGene("NCAM1")
res=getMutationStatsForGene("SHANK3")
res=getMutationStatsForGene("FOXO6")