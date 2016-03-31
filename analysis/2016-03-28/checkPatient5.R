library(synapseClient)
synapseLogin()

##patient 5 is missing the NF1 mutation 
high.tab<-read.table(synGet('syn5579055')@filePath,header=T,sep=' ',quote='"')
med.tab<-read.table(synGet('syn5586866')@filePath,header=T,sep=' ',quote='"')
low.tab<-read.table(synGet('syn5586903')@filePath,header=T,sep=' ',quote='"')

full.tab<-rbind(high.tab,med.tab,low.tab)
library(ggplot2)
full.tab$MinAllFreq=as.numeric(sapply(full.tab$GMAF,function(x) unlist(strsplit(as.character(x),split=':'))[2]))


ggplot(full.tab)+stat_bin(aes(MinAllFreq,fill=FILTER))

