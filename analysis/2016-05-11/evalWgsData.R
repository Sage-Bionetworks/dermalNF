##test wgs processing
require(parallel)
#lapply<<-function(x) mclapply(x,mc.cores=16)
mca<-function(x,FUN,...) mclapply(x,FUN,...,mc.cores=10)
assign('lapply',mca)
source("../../bin/WGSData_VarDict.R")
storeMutsForAllGenes(impact=c("HIGH"),0.05)
storeMutsForAllGenes(impact=c("MODERATE"),0.05)
storeMutsForAllGenes(impact=c("LOW"),0.05)
