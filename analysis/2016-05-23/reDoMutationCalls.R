###create updated list
source("../../dermalNF/bin/WGSData_VarDict.R")

res<-divideMAFFiles()
res.high<-divideMAFFiles(effect=c("HIGH"))