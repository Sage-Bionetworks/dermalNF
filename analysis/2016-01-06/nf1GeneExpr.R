##get NF1  expression per patient

source("../../bin/dermalNFData.R")

rna.seq<-rna_count_matrix(doNorm=TRUE)

rna.ann<-rna_annotations()

##now create data frame of NF1 gene expressiond data
library(ggplot2)

nf1.expr<-rna.seq['NF1',]
patients<-sapply(rna.ann[match(names(nf1.expr),rna.ann$synapseId),'patientId'],function(x) gsub('CT0+','',x))
tiss=rna.ann[match(names(nf1.expr),rna.ann$synapseId),'tissueId']

df<-data.frame(Expression=nf1.expr,Patient=patients,Tissue=tiss)
png('NF1_rnaExpr.png')
p<-ggplot(df)+geom_jitter(aes(Patient,Expression,colour=Patient),height=0)+ggtitle("NF1 Expression across patients")
print(p)
dev.off()

source("../../bin/WGSData.R")
##now can we somehow overlay? 
muts<-getMutationStatsForGene("NF1")

##get dnarnamatch
drmatch<-synTableQuery("SELECT Patient,DnaID,RnaID FROM syn5556216")@values

mutNum<-apply(df,1,function(x) {
  pats=which(as.character(muts$Patient)==as.character(x['Patient']))
  print(pats)
  if(length(pats)==0)
    return(-1)
  dids=intersect(which(drmatch$Patient==as.character(x[['Patient']])),which(drmatch$RnaID==as.numeric(as.character(x[['Tissue']]))))
  dd<-as.numeric(drmatch[dids,2])
  print(dd)
  if(is.na(dd))
    return(-1)
  else
    return(muts[intersect(pats,which(as.numeric(as.character(muts$Tissue))==dd)),'NumMutations'])
})
newdf=data.frame(df,SomaticMutations=as.factor(mutNum))

png('NF1_rnaExprWithMut.png')
p<-ggplot(newdf)+geom_jitter(aes(Patient,Expression,colour=Patient,shape=SomaticMutations),height=0)+ggtitle("NF1 Expression across patients")
print(p)
dev.off()

sdf<-subset(newdf,SomaticMutations!=-1)
png('NF1_rnaExprWithMut_matched.png')
p<-ggplot(sdf)+geom_jitter(aes(Patient,Expression,colour=Patient,shape=SomaticMutations),height=0)+ggtitle("NF1 Expression across patients")
print(p)
dev.off()

