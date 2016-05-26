##replot Concordances

library(data.table)
tab<-data.frame(fread('allSamples.tab',sep='\t',header=T,skip=58))
colnames(tab)<-c("CN","Discordance","NumberOfSites","AvMinDepth","Sample_i","Sample_j")

##now get annotations
require(synapseClient)
synapseLogin()
vcfs=synapseQuery("select name,id,patientID,tissueID from entity where parentId=='syn5522788'")
vcfs$Patient=sapply(vcfs$entity.patientID,function(x) gsub("CT0+",'',x))
vcfs$Base=sapply(vcfs$entity.name,function(x) gsub('.vcf','',x,fixed=T))
vcfs$Tissue=sapply(vcfs$entity.tissueID,function(x) gsub("[0-9]+","Tumor",x))
tab$Patient_i=vcfs$Patient[match(tab$Sample_i,vcfs$Base)]
tab$Patient_j=vcfs$Patient[match(tab$Sample_j,vcfs$Base)]
tab$Tissue_i=vcfs$Tissue[match(tab$Sample_i,vcfs$Base)]
tab$Tissue_j=vcfs$Tissue[match(tab$Sample_j,vcfs$Base)]

require(reshape2)
require(pheatmap)

atab<-tab
atab$Sample_i=tab$Sample_j
atab$Sample_j=tab$Sample_i
atab$Patient_i=tab$Patient_j
atab$Patient_j=tab$Patient_i
atab$Tissue_i=tab$Tissue_j
atab$Tissue_j=tab$Tissue_i
newtab=rbind(atab,tab)

dmat<-acast(newtab,Sample_j~Sample_i,value.var="Discordance")
itab<-unique(newtab[,c(5,7,9)])
ivals=data.frame(Tissue=itab$Tissue_i,Patient=itab$Patient_i)
rownames(ivals)<-itab$Sample_i

jtab<-unique(newtab[,c(6,8,10)])

jvals=data.frame(Tissue=jtab$Tissue_j,Patient=jtab$Patient_j)
rownames(jvals)<-jtab$Sample_j

pheatmap(dmat,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampleDiscordance.png')
pheatmap(dmat,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampleDiscordance.pdf')

numsites=acast(newtab,Sample_j~Sample_i,value.var='NumberOfSites')
pheatmap(numsites,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampNumOfsites.png')
pheatmap(numsites,annotation_row=ivals,cellheight=10,cellwidth=10,annotation_col=jvals,file='sampNumOfsites.pdf')


for (file in c('sampleDiscordance','sampNumOfsites')){
  for(suff in c('.pdf','.png')){  
    fname=paste(file,suff,sep='')
    synStore(File(fname,parentId='syn5669832'),
                 used=list(list(entity='syn5669991',wasExecuted=FALSE),
                           list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-09/replotConcordances.R',wasExecuted=TRUE)))
  }
}

