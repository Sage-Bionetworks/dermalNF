##begin to work on figures

source("../../bin/dermalNFData.R")

##figure 4 is the easiest, starting wtiht hat

library(ggplot2)
#library(dplyr)

##figure 3 focuses on SNP data
cnv_annotes=cnv_annotations()
library(parallel)
all.cnv=cnv_unprocessed()
lrr<- do.call("rbind", mclapply(all.cnv, function(i) as.data.frame(i[,2:3]),mc.cores=4))
baf<- do.call("rbind", mclapply(all.cnv, function(i) as.data.frame(i[,c(2,4)]),mc.cores=4))

cnv_annotes$patient=sapply(as.character(cnv_annotes$patientId),function(x) gsub("CT0+","",x))
##add in patient identifier
lfnames=sapply(lrr$Sample.ID,paste,'Final.csv',sep='_')
lrr$Patient=cnv_annotes$patient[match(lfnames,cnv_annotes$File)]
lrr$Samples=paste('Patient',lrr$Patient,'Tissue',cnv_annotes$tissueID[match(lfnames,cnv_annotes$File)],sep='_')
lrr$tissueType=cnv_annotes$tissueType[match(lfnames,cnv_annotes$File)]

fnames=sapply(baf$Sample.ID,paste,'Final.csv',sep='_')
baf$Patient=cnv_annotes$patient[match(fnames,cnv_annotes$File)]
baf$tissueType=cnv_annotes$tissueType[match(fnames,cnv_annotes$File)]
baf$Samples=paste('Patient',baf$Patient,'Tissue',cnv_annotes$tissueID[match(lfnames,cnv_annotes$File)],sep='_')


##figure 5 - plot of proteomics values?
pl=ggplot(lrr,aes(y=Log.R.Ratio,x=tissueType))+geom_violin(aes(fill=Patient,colour=tissueType)) + coord_flip()
pb=ggplot(lrr,aes(y=Log.R.Ratio,x=tissueType))+geom_boxplot(aes(fill=Patient,colour=tissueType)) + coord_flip()

pl2=ggplot(lrr,aes(y=Log.R.Ratio,x=Samples))+geom_violin(aes(fill=Patient,colour=tissueType)) + coord_flip()
pb2=ggplot(lrr,aes(y=Log.R.Ratio,x=Samples))+geom_boxplot(aes(fill=Patient,colour=tissueType)) + coord_flip()

pdf('rotated_LrrPlot.pdf')#,height=800)

print(pl)
print(pb)
print(pl2)
print(pb2)
dev.off()

pb=ggplot(baf,aes(y=B.Allele.Freq,x=tissueType))+geom_violin(aes(fill=Patient,colour=tissueType))+coord_flip()
pb2=ggplot(baf,aes(y=B.Allele.Freq,x=Samples))+geom_violin(aes(fill=Patient,colour=tissueType))+coord_flip()

pdf('rotated_violinBafPlot.pdf')#,height=800)
print(pb)
print(pb2)
dev.off()

snpqc='syn5669811'
scripturl='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-07/figure_3_draft.R'

for(file in c('rotated_violinBafPlot.pdf','rotated_LrrPlot.pdf')){
  synStore(File(file,parentId=snpqc),
           used=list(list(url=scripturl,wasExecuted=TRUE)))
}
