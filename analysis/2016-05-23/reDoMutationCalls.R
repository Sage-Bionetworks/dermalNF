###create updated list
source("../../bin/WGSData_VarDict.R")

<<<<<<< HEAD
#res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'))  #)divideMAFfiles()
#res.high<-storeMutsForAllGenes(impact=c("HIGH"))

res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'),pvalthresh=0.1)  #)divideMAFfiles()
res.high<-storeMutsForAllGenes(impact=c("HIGH"),pvalthresh=0.1)
=======
res<-storeMutsForAllGenes(impact=c('HIGH','MODERATE','LOW'),0.1)  #)divideMAFfiles()
#res.high<-storeMutsForAllGenes(impact=c("HIGH"))
>>>>>>> e494c8e667968edb9332fefedd886e7bc6aaff6e
