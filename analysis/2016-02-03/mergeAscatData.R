###take ASCAT data and merge into matrices

ascat.dir='../../ASCAT/'
ascat.dir='~tyu/projects/DERMALNF/'
source("../../bin/dermalNFData.R")
##get teh CNV annotations
annotes=cnv_annotations()

##get all files in directory
allfiles=list.files(ascat.dir)
logrfiles=allfiles[grep('LogR',allfiles)]
baffiles=allfiles[grep('BAF',allfiles)]

logr.mat=mclapply(logrfiles,function(x) read.table(paste(ascat.dir,x,sep=''),sep='\t',row.names=1))

sids=sapply(logrfiles,function(x){
    fp=gsub('.txt','',unlist(strsplit(x,split='_X'))[2])
    idx=grep(gsub('.','-',fp,fixed=T),annotes$File)
    sid=paste(annotes$synapseId[idx],annotes$patientId[idx],annotes$tissueId[idx],sep='_')
    sid})
logr.mat<-as.data.frame(logr.mat)
colnames(logr.mat)<-sids

write.table(logr.mat,'ascat_tumor_normal_logR.txt',quote=F,row.names=T,col.names=T)

baf.mat=mclapply(baffiles,function(x) {
    tab=read.table(paste(ascat.dir,x,sep=''),sep='\t')
    ntab=cbind(tab,rep(x,nrow(tab)))
    colnames(ntab)=c("SNP","BAF","FILE")
    ntab})

require(reshape2)
require(plyr)
baf.df=ldply(baf.mat,rbind)
                                        #now reshape as matrix
baf.sq=acast(baf.df,SNP~FILE,value.var='BAF',fill=NA)


sids=sapply(baffiles,function(x){
    fp=gsub('.txt','',unlist(strsplit(x,split='_X'))[2])
    idx=grep(gsub('.','-',fp,fixed=T),annotes$File)
    sid=paste(annotes$synapseId[idx],annotes$patientId[idx],annotes$tissueId[idx],sep='_')
    sid})

colnames(baf.sq)<-sids[colnames(baf.sq)]

write.table(baf.sq,'ascat_tumor_normal_BAF.txt',quote=F,row.names=T,col.names=T)

require(synapseClient)
lf=File('ascat_tumor_normal_logR.txt',parentId='syn5049702')
synStore(lf,executed=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-03/mergeAscatData.R')))

af=File('ascat_tumor_normal_BAF.txt',parentId='syn5049702')

synStore(af,executed=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-03/mergeAscatData.R')))
