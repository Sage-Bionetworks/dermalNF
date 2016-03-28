library(synapseClient)
synapseLogin()
high.tab<-read.table(synGet('syn5578970')@filePath,header=T,sep=' ',quote='"')
med.tab<-read.table(synGet('syn5586790')@filePath,header=T,sep=' ',quote='"')
low.tab<-read.table(synGet('syn5586829')@filePath,header=T,sep=' ',quote='"')

full.tab<-rbind(high.tab,med.tab,low.tab)
library(ggplot2)
full.tab$MinAllFreq=as.numeric(sapply(full.tab$GMAF,function(x) unlist(strsplit(as.character(x),split=':'))[2]))


ggplot(full.tab)+stat_bin(aes(MinAllFreq,fill=FILTER))

