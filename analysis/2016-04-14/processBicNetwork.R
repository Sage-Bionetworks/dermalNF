library(synapseClient)
synapseLogin()

load(synGet('syn5908832')@filePath)

net<-read.csv('dermalNFrankconsensus.csv')
ppi<-cbind(as.character(net$var1),rep('pp',nrow(net)),as.character(net$var2))
write.table(ppi,'dermalNFrankConsensus.sif',sep='\t',quote=F,row.names=F,col.names=F)
edges<-which(bicNetworks$sparrow$network,arr.ind=T)

sp.ppi<-cbind(rownames(bicNetworks$sparrow$network)[edges[,1]],rep('pp',nrow(edges)),rownames(bicNetworks$sparrow$network)[edges[,2]])
write.table(sp.ppi,'dermalNFsparrowNetwork.sif',sep='\t',quote=F,row.names=F,col.names=F)

mirsAndSnosAndLincs<-c(grep('MIR',sp.ppi[,1]),grep('MIR',sp.ppi[,3]),
                       grep('SNO',sp.ppi[,1]),grep('SNO',sp.ppi[,3]),
                       grep('LINC',sp.ppi[,1]),grep('LINC',sp.ppi[,3]))
write.table(sp.ppi[-unique(mirsAndSnosAndLincs),],'dermalNFsparrowNetworkNosmRNA.sif',sep='\t',quote=F,row.names=F,col.names=F)