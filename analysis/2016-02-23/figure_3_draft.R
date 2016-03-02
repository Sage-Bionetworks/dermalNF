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
lrr$tissueType=cnv_annotes$tissueType[match(lfnames,cnv_annotes$File)]

fnames=sapply(baf$Sample.ID,paste,'Final.csv',sep='_')
baf$Patient=cnv_annotes$patient[match(fnames,cnv_annotes$File)]
baf$tissueType=cnv_annotes$tissueType[match(fnames,cnv_annotes$File)]

##figure 5 - plot of proteomics values? 
pl=ggplot(lrr,aes(y=Log.R.Ratio,x=Sample.ID))+geom_violin(aes(fill=Patient,colour=tissueType))+coord_flip()
png('rotated_violinLrrPlot.png',height=800)

print(pl)
dev.off()

pb=ggplot(baf,aes(y=B.Allele.Freq,x=Sample.ID))+geom_violin(aes(fill=Patient,colour=tissueType))+coord_flip()
png('rotated_violinBafPlot.png',height=800)
print(pb)
dev.off()

snpqc='syn5669811'
scripturl='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-23/figure_3_draft.R'

for(file in c('rotated_violinBafPlot.png','rotated_violinLrrPlot.png')){
  synStore(File(file,parentId=snpqc),
           used=list(
                     list(url=scripturl,wasExecuted=TRUE)))
}

