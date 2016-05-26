#doing more exploration of RNA-Seq data. Starting with recapitulating analysis on 01/06
source("../../bin/dermalNFData.R")

##get mutational data too
source("../../bin/WGSData.R")
##now can we somehow overlay? 
muts<-getMutationStatsForGene("NF1")

gene.cor=''
samp.cor=''


##first plot histogram of correlations


##now get cpount files, then plot those for NF1 with sample data
rna.seq<-rna_count_matrix(doNorm=TRUE,minCount=2)
rna.fpkm<-rna_fpkm_matrix()
rna.ann<-rna_annotations()

##now create data frame of NF1 gene expressiond data
library(ggplot2)

nf1.count.expr<-rna.seq['NF1',]
nf1.fpkm.expr<-rna.fpkm['NF1',]
##Not enough samples for FPKM- redo once. we have 
patients<-sapply(rna.ann[match(names(nf1.count.expr),rna.ann$synapseId),'patientId'],function(x) gsub('CT0+','',x))
tiss=rna.ann[match(names(nf1.count.expr),rna.ann$synapseId),'tissueId']
counts.df<-data.frame(Expression=nf1.count.expr,Patient=as.character(patients),Tissue=as.numeric(as.character(tiss)))


####what other analysis can we do? 
#let's revisit the NF1 expression/snp analysis

drmatch<-synTableQuery("SELECT Patient,DnaID,RnaID FROM syn5556216")@values

mutNum<-apply(counts.df,1,function(x) {
  ##which mutations belong to patient?
  pats=match(as.character(x[["Patient"]]),as.character(muts$Patient))
  dids=intersect(which(drmatch$RnaID==x[['Tissue']]),which(as.character(drmatch$Patient)==as.character(x[['Patient']])))
  dd<-as.numeric(drmatch[dids,2])
  ##
  if(length(pats)==0)
    return(c(Patient=as.character(x[['Patient']]),Somatic=-1,Germline=-1))
  if(is.na(dd))
    return(c(Patient=as.character(x[['Patient']]),Somatic=-1,Germline=-1))
  else{
    am=muts[intersect(pats,which(as.numeric(as.character(muts$Tissue))==dd)),]
    return(c(Patient=as.character(x[['Patient']]),Somatic=length(which(am$MutationType=="Somatic")),Germline=length(which(am$MutationType=='Germline'))))

  }
})
newdf=data.frame(counts.df,SomaticMutations=as.factor(mutNum['Somatic',]))
newdf$GermlineMutations=as.factor(sapply(counts.df$Patient,function(x)length(unique(subset(subset(muts,Patient==as.character(x)),MutationType=='Germline')$Position))))

png('NF1_somatic_germline_muts_expresion.png')
p<-ggplot(newdf)+geom_jitter(aes(Patient,Expression,colour=Patient,shape=SomaticMutations,size=GermlineMutations),height=0)+ggtitle("NF1 Expression across patients")
print(p)
dev.off()

png('NF1_somatic_germline_rank_muts_expresion.png')
p<-ggplot(newdf)+geom_point(aes(rank(Expression),Expression,colour=Patient,shape=SomaticMutations,size=GermlineMutations))+ggtitle("NF1 Expression across patients")
print(p)
dev.off()

##let's try to build linear model?
allcoeffs<-t(apply(rna.seq,1,function(x){
  adf=newdf
  adf$Expression=x
  coefficients(summary(lm(Expression~GermlineMutations+SomaticMutations+GermlineMutations*SomaticMutations,adf)))[-1,4]
}))
cor.coeffs<-apply(allcoeffs,2,p.adjust)

##what genes correlat with NF1 mutation?
mgenes=unique(rownames(which(cor.coeffs<0.05,arr.ind=T)))

if(is.null(mgenes))
  mgenes=unique(rownames(which(allcoeffs<0.005,arr.ind=T)))

sapply(mgenes,function(x){
  adf<-newdf
  adf$Expression=rna.seq[x,]
  png(paste(x,'_somatic_germline_mixed_model_ranked_expresion.png',sep=''))
  p<-ggplot(adf)+geom_point(aes(rank(Expression),Expression,colour=Patient,shape=SomaticMutations,size=GermlineMutations))+ggtitle(paste(x,"Expression across patients"))
  print(p)
  dev.off()
})

#now look only at germline
germ.coeffs=t(apply(rna.seq,1,function(x){
  adf=newdf
  adf$Expression=x
  coefficients(summary(lm(Expression~GermlineMutations,adf)))[-1,4]
}))
cor.germ<-p.adjust(germ.coeffs[1,])
#now look at somatic
som.coeffs=t(apply(rna.seq,1,function(x){
  adf=newdf
  adf$Expression=x
  adf<-subset(adf,SomaticMutations!=-1)
  coefficients(summary(lm(Expression~SomaticMutations,adf)))[-1,4]
}))
cor.som<-p.adjust(som.coeffs[1,])

anymut.coffs=t(apply(rna.seq,1,function(x){
  adf=newdf
  adf$Expression=x
  adf<-subset(adf,SomaticMutations!=-1)
  adf$EitherMutation=as.numeric(as.character(adf$SomaticMutations))+as.numeric(as.character(adf$GermlineMutations))
  coefficients(summary(lm(Expression~EitherMutation,adf)))[-1,4]
}))
any.corr<-p.adjust(anymut.coffs[1,])
#again, nothing passes.

mgenes=rownames(rna.seq)[order(anymut.coffs)[1:30]]
sapply(mgenes,function(x){
  adf<-newdf
  adf$Expression=rna.seq[x,]
  png(paste(x,'_somatic_germline_combined_model_ranked_expresion.png',sep=''))
  p<-ggplot(adf)+geom_point(aes(rank(Expression),Expression,colour=Patient,shape=SomaticMutations,size=GermlineMutations))+ggtitle(paste(x,"Expression across patients"))
  print(p)
  dev.off()
})

##no real correlation between various trancscripts and nF1 mutation status

nf1.cors=apply(rna.seq,1,function(x)
  cor(x,rna.seq['NF1',]))

msi2.cors=apply(rna.seq,1,function(x)
  cor(x,rna.seq['MSI2',]))


muts<-getMutationStatsForGene("SEC31B")

muts<-getMutationStatsForGene("")
topstats=sapply(names(sort(nf1.cors,decreasing=T)[1:30]),getMutationStatsForGene)
nextstats=sapply(names(sort(nf1.cors,decreasing=T)[31:50]),getMutationStatsForGene)