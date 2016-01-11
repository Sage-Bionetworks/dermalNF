##let's measure mRNA variation across/within samples

#norm counts
source('../../bin/dermalNFData.R')

counts<-rna_count_matrix(doNorm=TRUE, minCount=2)

rna_annotes<-rna_annotations()

##now create vector of patient factors

Patients<-sapply(rna_annotes$patientId,function(x) gsub('CT0+','',x))
names(Patients)<-rna_annotes$synapseId
df<-data.frame(Patients=as.factor(Patients),t(counts[,match(rna_annotes$synapseId,colnames(counts))]))

genenames<-colnames(df)[-1]

#do linear model of each gene
all.lm<-lapply(genenames,function(x){
  res<-summary(lm(as.formula(paste(x,'Patients',sep='~')),data=df))
  return(res)
})

res.coeff<-sapply(all.lm,function(x) coefficients(x)[,4])
res.ps<-sapply(all.lm,function(x){
  fs=x$fstatistic
  pf(fs[1],fs[2],fs[3],lower.tail=F)
})

adj=p.adjust(res.ps,'BH')

##now assess results

sig.genes<-genenames[which(adj<0.1)]
fname=paste('geneExpressionOf',length(sig.genes),'patSpecGenes.png',sep='')
pheatmap(log2(counts[sig.genes,]),cellwidth=10,cellheight=10,annotation_col=data.frame(Patient=Patients),filename=fname)

write(sig.genes,file='genesAssociatedWithPatients.txt')
write(rownames(counts),file='bgGenes.txt')
write(genenames[order(res.ps,decreasing =F)],file='modelAssociatedGenesRanked.txt')