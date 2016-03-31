####compare TCGA expression data with dermal data

source('../../bin/TcgaExpressionData.R')
source('../../bin/TcgaMutationalData.R')

##now get dermalNF data

load('combinedExprPC.rda')
library(ggbiplot)
#ggbiplot(pc,var.axes=F)

computeCentroidsForDisease<-function(pc,patlist){
    
  
}

