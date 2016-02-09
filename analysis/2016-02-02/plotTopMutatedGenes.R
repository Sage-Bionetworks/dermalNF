###goal is to plot the most mutated genes across all samples

som.muts<-read.table('../2016-01-27/somaticCancerGeneMutationsInDermals.tsv',sep='\t',header=T)
germ.muts<-read.table('../2016-01-27/germlineCancerGeneMutationsInDermals.tsv',sep='\t',header=T)


##we forgot to ascertain copy number loss! 
#let's do that now and put into a heatmap to characterize
source("../../bin/dermalNFData.R")
cancer.genes=read.table(synGet('syn5610781')@filePath,sep=',',header=T)
cnv_seg=cnv_segmented_by_gene()
##now filter that for cancer genes
idx=setdiff(match(as.character(cancer.genes[,1]),rownames(cnv_seg)),NA)
cancer.cnv<-cnv_seg[idx,]
library(pheatmap)

sampDat<-data.frame(Tissue=tissueType,Patient=patients)
pheatmap(cancer.cnv,cellheight=10,cellwidth=10,annotation_col=sampDat,
         file='copyNumberLogRCancerGenes.png')

##now restrict for some greater loss
lower.vals=apply(cancer.cnv,1,function(x) any(x<(-0.1)))

pheatmap(cancer.cnv[which(lower.vals),],cellheight=10,cellwidth=10,
         annotation_col=sampDat,filename='copyNumberLogRatioBelow0.1.png')

##now can we roughly compute somatic copy number? 
fulldat<-sampDat
fulldat$synId=rownames(sampDat)

fulldat=subset(fulldat,!Patient%in%c(13,8))
cndiff=apply(subset(fulldat,Tissue=='Tissue'),1,function(x){
  print(x)
  tiss=cancer.cnv[,x[['synId']]]
  blood=cancer.cnv[,subset(subset(fulldat,Tissue=='Blood'),Patient==x[['Patient']])[['synId']]]
  res=tiss-blood
  names(res)=rownames(cancer.cnv)
  res
})

##now can we plot
pheatmap(cndiff,cellheight=10,cellwidth=10,annotation_col=sampDat,
         file='somaticCopyNumberLogRCancerGenes.png')

##now restrict for some greater loss
lower.vals=apply(cndiff,1,function(x) any(x<(-0.1)))

pheatmap(cndiff[which(lower.vals),],cellheight=10,cellwidth=10,
         annotation_col=sampDat,filename='somaticCopyNumberLogRatioBelow0.1.png')

##now create visualization of WGS data!!!
require(dplyr)
require(reshape2)
som.muts=subset(som.muts,Mutation_Type!='Silent')
germ.muts=subset(germ.muts,Mutation_Type!='Silent')

som.by.pat=group_by(som.muts,Hugo_Symbol) %>% summarize(n_distinct(Sample_ID))
germ.by.pat=group_by(germ.muts,Hugo_Symbol) %>% summarize(n_distinct(Sample_ID))

mut.mat=acast(som.muts,Hugo_Symbol~Sample_ID,value.var = "Mutation_Type",fun.aggregate = length)

mut.mat.no.10<-mut.mat[,-grep('patient_10',colnames(mut.mat))]
nz=which(apply(mut.mat.no.10,1,function(x) any(x>0)))
pheatmap(mut.mat.no.10[nz,])


mut.mat.no.8.10<-mut.mat.no.10[,-grep('patient_8',colnames(mut.mat.no.10))]
nz=which(apply(mut.mat.no.8.10,1,function(x) any(x>0)))
pheatmap(mut.mat.no.8.10[nz,])

com.muts=which(apply(mut.mat,1,function(x) length(which(x!=0)))>5)
pheatmap(log10(1+mut.mat[com.muts,]))

gmut.mat=acast(germ.muts,Hugo_Symbol~Sample_ID,value.var = "Mutation_Type",fun.aggregate = length)
com.muts=which(apply(gmut.mat,1,function(x) length(which(x!=0)))>8)
pheatmap(log10(1+gmut.mat[com.muts,]))
