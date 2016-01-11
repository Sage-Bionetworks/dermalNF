##apply gene expression profiles to GSEA/MSIGDB using GSVA
library(GSVA)

source('../../bin/dermalNFData.R')


library(data.table)
library(GSEABase)
library(GSVAdata)
#cut and paste from vignette
#cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,
#                            +                            min.sz=10, max.sz=500, verbose=TRUE)$es.obs,
#        +                            dir=cacheDir, prefix=cachePrefix)

data(c2BroadSets)

#> gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)$es.obs
rna.data<-rna_count_matrix()
#get entrez ids for RNA count matrix
eids<-as.data.frame(fread('../../data/HugoGIDsToEntrez_DAVID.txt',sep='\t',header=T))

eid.match=eids[match(rownames(rna.data),eids[,1]),2]
eid.rows=which(!is.na(eid.match))
ent.rna<-rna.data[eid.rows,]
rownames(ent.rna)=eid.match[eid.rows]

es<-gsva(ent.rna,c2BroadSets)

##now get patient annotaions
pats=rna.annot$patientId
tiss=rna.annot$tissueID

newnames=paste(pats,tiss)
names(newnames)<-rna.annot$synapseId

#let's cluster, see how that goes? 
res=es$es.obs
colnames(res)<-newnames[colnames(res)]
plot(hclust(dist(t(res)),method='ward.D2'))

res<-es$es.obs
vars<-apply(res,1,var,na.rm=T)
mostvar<-res[order(vars,decreasing=T)[1:100],]
names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10,annotation_col=data.frame(Patient=pats) ,filename='mostVariablePathways_GSVA.png')

