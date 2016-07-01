##now that we have ascat data processed, we can 


source("../../bin/dermalNFData.R")


ascat.vals<-ascat_segments(recalc=TRUE,annot=NA,byval='gene',metric='median')

scripts=c('https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-06-30/buildAscatAnalysis.R','https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/dermalNFData.R')
synvals='syn6182422'


logr<-ascat.vals$LRR
baf<-ascat.vals$BAF
logr.seg<-ascat.vals$LRR.seg
baf.seg<-ascat.vals$BAF.seg

write.table(logr,file='ASCAT_calculated_smoothedLRRGeneVals.tsv',sep='\t',col.names=T,row.names=F,quote=F)
write.table(baf,file='ASCAT_calculated_smoothedBAFGeneVals.tsv',sep='\t',col.names=T,row.names=F,quote=F)
write.table(logr.seg,file='ASCAT_calculated_smoothedLRRVals.seg',sep='\t',col.names=T,row.names=F,quote=F)
write.table(baf.seg,file='ASCAT_calculated_smoothedBAFVals.seg',sep='\t',col.names=T,row.names=F,quote=F)

synStore(File('ASCAT_calculated_smoothedLRRVals.seg',parentId='syn6182411'),executed=list(list(url=scripts[1]),list(url=scripts[2])),
         used=list(list(entity=synvals)))

synStore(File('ASCAT_calculated_smoothedBAFVals.seg',parentId='syn6182411'),executed=list(list(url=scripts[1]),list(url=scripts[2])),
         used=list(list(entity=synvals)))

synStore(File('ASCAT_calculated_smoothedLRRGeneVals.tsv',parentId='syn6182411'),executed=list(list(url=scripts[1]),list(url=scripts[2])),
         used=list(list(entity=synvals)))

synStore(File('ASCAT_calculated_smoothedBAFGeneVals.tsv',parentId='syn6182411'),executed=list(list(url=scripts[1]),list(url=scripts[2])),
         used=list(list(entity=synvals)))