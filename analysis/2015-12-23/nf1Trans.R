##plot NF1 transcripts across dermalNF samples
source("../../bin/dermalNFData.R")

rna.counts<-rna_count_matrix(doNorm=TRUE)

nf<-which(rownames(rna.counts)=='NF1')
annotes<-rna_bam_annotations()
annotes=annotes[which(!is.na(annotes$patientID)),]
samps<-apply(annotes,1,function(x) paste(x[c(1,4)],collapse='_samp_'))
samps<-sapply(samps,function(x)(gsub("CT0+","Patient",x)))
colnames(rna.counts)<-samps

library(ggplot)
df<-data.frame(Counts=rna.counts[nf,],Patient=sapply(samps,function(x) unlist(strsplit(x,split='_'))[1]))
#plot(rna.counts[nf,])
png("normalizedNF1Counts.png")
p<-ggplot(df)+geom_boxplot(aes(Patient,Counts,border=Patient))+geom_jitter(aes(Patient,Counts,color=Patient))
print(p)
dev.off()