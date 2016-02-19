##goal of this script is to visualize mutational data (possibly alongside copy number?)


library(synapseClient)

all.muts<-read.table(synGet('syn5611520')@filePath,sep='\t',header=T)

collectGermlineMuts<-function(){
  
  
}

collectSomaticMuts<-function(){
  
  
}