###create updated list
source("../../bin/WGSData_VarDict.R")

#res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'))  #)divideMAFfiles()
#res.high<-storeMutsForAllGenes(impact=c("HIGH"))

res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'),pvalthresh=0.1)  #)divideMAFfiles()
res.high<-storeMutsForAllGenes(impact=c("HIGH"),pvalthresh=0.1)
