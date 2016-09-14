##re-quantify proteomic data using proteome discover before proteoIQ
source("../../bin/dermalNFData.R")

processProtFile<-function(synid){
  tab<-read.table(synGet(synid)@filePath,sep='\t',header=T,as.is=T)
  gene.quants<-apply(tab,1,function(x){
      gene=gsub('gene:','',unlist(strsplit(x[10],split=' '))[4])
      
  })
}