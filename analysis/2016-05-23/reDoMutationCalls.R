###create updated list
source("../../bin/WGSData_VarDict.R")

res<-storeMutsForAllGenes(effect=c('HIGH','MODERATE','LOW')  #)divideMAFfiles()
res.high<-storeMutsForAllGenes(effect=c("HIGH"))
