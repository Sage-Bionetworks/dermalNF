###collect copy number alterations, map to getnes
source("../../bin/dermalNFData.R")

##get mapping file
lrr=read.table(synGet('syn5622705')@filePath)

baf=read.table(synGet('syn5622683')@filePath)

annot <- snp_annotation_data()
annot.idx=match(rownames(lrr),annot$Name)

all.pos<-annot$MapInfo[annot.idx]
all.chr<-annot$Chr[annot.idx]
is.autosome <- as.character(all.chr) %in% as.character(1:22)

rm(annot)


library(DNAcopy)
cna <- CNA(lrr[is.autosome,], as.character(all.chr)[is.autosome], all.pos[is.autosome], data.type="logratio",colnames(lrr))


smoothed.cna <- smooth.CNA(cna)
segment.smoothed.cna <- segment(smoothed.cna, verbose=1)

segment.smoothed.cna.sundo <- segment(smoothed.cna, undo.splits="sdundo",undo.SD=2,verbose=1)

write.table(segment.smoothed.cna$output, file="dermal_nf_SOMATIC_cbs_noundo.seg",
            sep="\t",quote=FALSE,row.names=F)
write.table(segment.smoothed.cna.sundo$output, file="dermal_nf_SOMATIC_cbs_undosd2.seg",
            sep="\t",quote=FALSE,row.names=F)


##now we can do the gene mapping as well.
if(!exists("geneInfo"))
  geneInfo<-read.table('../../data/hg19_geneInfo.txt')

segdat<-segment.smoothed.cna.sundo$output

library(CNTools)
cnseg <- CNSeg(segdat)
metric='median'
byval='gene'

rdseg <- getRS(cnseg, by = byval,geneMap=geneInfo, imput = FALSE, XY = FALSE, what =metric)
segM <- rs(rdseg)
#nzvals<-which(apply(segM,1,function(x) any(as.numeric(x[-c(1:5)])<thresh)))
#nzM<-segM[nzvals,]

otherSeg=cnv_segmented(TRUE)
cnseg2 <- CNSeg(otherSeg)
rdseg2 <- getRS(cnseg2, by = byval,geneMap=geneInfo, imput = FALSE, XY = FALSE, what =metric)
segM2 <- rs(rdseg2)

g1=which(segM$genename=='NF1')
g2=which(segM2$genename=="NF1")

sort(c(colSums(segM2[g2,-c(1:5)]),colSums(segM[g1,-c(1:5)])))
     

