###create updated list
source("../../bin/WGSData_VarDict.R")

res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'))  #)divideMAFfiles()
res.high<-storeMutsForAllGenes(impact=c("HIGH"))
