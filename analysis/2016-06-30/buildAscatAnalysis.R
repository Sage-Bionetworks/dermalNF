##now that we have ascat data processed, we can 


source("../../bin/dermalNFData.R")


ascat.vals<-ascat_segments(annot=NA,byval='gene',metric='median')

scripts=c('','')
synvals='syn6182422'