###collect copy number alterations, map to getnes
source("../../bin/dermalNFData.R")

##get mapping file
annot <- snp_annotation_data()
is.autosome <- as.character(annot$Chr) %in% as.character(1:22)

all.pos<-annot$MapInfo
all.chr<-annot$Chr

rm(annot)

lrr=read.table(synGet('')@filePath)

cna <- CNA(lrr[is.autosome,], as.character(all.chr)[is.autosome], all.pos[is.autosome], data.type="logratio",names(sample.data))


smoothed.cna <- smooth.CNA(cna)
segment.smoothed.cna <- segment(smoothed.cna, verbose=1)

segment.smoothed.cna.sundo <- segment(smoothed.cna, undo.splits="sdundo",undo.SD=2,verbose=1)

write.table(segment.smoothed.cna$output, file="dermal_nf_SOMATIC_cbs_noundo.seg",
            sep="\t",quote=FALSE,row.names=F)
write.table(segment.smoothed.cna.sundo$output, file="dermal_nf_SOMATIC_cbs_undosd2.seg",
            sep="\t",quote=FALSE,row.names=F)
