####compare TCGA expression data with dermal data

source('../../bin/TcgaExpressionData.R')
source('../../bin/TcgaMutationalData.R')

##now get dermalNF data

source("../../bin/dermalNFData.R")
rna.counts<-rna_count_matrix(doLogNorm=T,minCount=2)
#rna.fpkm<-rna_fpkm_matrix()

expr.pats<-toPatientId(colnames(alldat))
gene.symbs<-sapply(alldat[,1],function(x) unlist(strsplit(x,split='|',fixed=T))[1])
fpkm.idx <- match(rownames(rna.fpkm),gene.symbs)
#counts.idx <-match(rownames(rna.counts),gene.symbs)

disdat<-alldat[,intersect(colnames(alldat),unlist(tumsByDis))]
##can we just get all cors
full.mat<-cbind(rna.counts[!is.na(counts.idx),],disdat[counts.idx[!is.na(counts.idx)],])

full.cor<-cor(full.mat,use='pairwise.complete.obs',method='spearman')


write.table(full.cor,file='panTcgaDermalCountsSpearmanCors.txt',sep='\t',row.names=T,col.names=T,quote=F)
synStore(File('panTcgaDermalCountsSpearmanCors.txt',parentId='syn5821631'),
         used=list(list(entity='syn4311114'),list(entity='syn3281840'),
                   list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-23/panExpressionComparison.R',wasExecuted=TRUE)))

##now compute the same for fpkm
full.mat<-cbind(rna.fpkm[!is.na(fpkm.idx),],disdat[fpkm.idx[!is.na(fpkm.idx)],])

full.cor<-cor(full.mat,use='pairwise.complete.obs',method='spearman')


write.table(full.cor,file='panTcgaDermalFPKMSpearmanCors.txt',sep='\t',row.names=T,col.names=T,quote=F)
synStore(File('panTcgaDermalFPKMSpearmanCors.txt',parentId='syn5821631'),
         used=list(list(entity='syn4311114'),list(entity='syn3281840'),
                   list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-23/panExpressionComparison.R',wasExecuted=TRUE)))

