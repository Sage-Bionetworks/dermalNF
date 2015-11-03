###
### DermalNF CVN analysis
### Goal of this script is to collect the CNV data from the OMNI SNP platform and determine
### 1- that the data is sound, the log ratios are distributed as expected
### 2- if there is loss in the NF1 region

library(data.table)
library(DNAcopy)
library(CNTools)
library(synapseClient)

#first login to synapse
synapseLogin()

                                        #read in the annotation file. this can be big.
anndata<-synGet('syn5005069')
annot <- as.data.frame(fread(anndata@filePath,sep=",",header=T))

##get regions that are within our general region: chr17:29000019 to 30427403
chr17.snps=annot[grep('chr17',annot$chrpos),]

all.pos<-as.numeric(sapply(annot$chrpos,function(x) unlist(strsplit(x,split='_'))[2]))
all.chr<-sapply(annot$chrpos,function(x) unlist(strsplit(x,split='_'))[1])

pos<-all.pos[grep('chr17',annot$chrpos)]
chr17.snps=chr17.snps[intersect(which(pos>29000019),which(pos<30427403)),]


#now collect samples from synapse
##all SNP files, need to add metadata...
snpfiles=synapseQuery('SELECT id,name,Patient_ID,Tissue_Type,Tissue_ID FROM entity where parentId=="syn5004874"')

snpfiles<-snpfiles[grep("Final.csv",snpfiles[,1]),]
sample.names<-sapply(snpfiles[,1],function(x) gsub('_Final.csv','',unlist(strsplit(x,split='-'))[3]))


sample.data<-lapply(snpfiles[,2],function(synid){
    fname=synGet(synid)
    data <- as.data.frame(fread(fname@filePath,sep=",",header=T))
    ad<-data[match(annot$Name,data$'SNP.Name'),]
    return(ad)
})

#remove x/y here
is.autosome <- as.character(annot$chr) %in% as.character(1:22)
#map lrr and bafs for each
lrr <- do.call("cbind", lapply(sample.data, function(x) x$"Log.R.Ratio"))
baf <- do.call("cbind", lapply(sample.data, function(x) x$"B.Allele.Freq"))
##add variable to signify if SNP is in region of interest
in.region<-rep(FALSE,nrow(lrr))

in.region[which(annot$chrpos%in%chr17.snps$chrpos)]<-TRUE


#how can we plot baf only?
library(ggplot2)
##create ggplot-associated data.frame
list.of.lists<-c(LogRatio=c(),BAlleleFreq=c(),Sample=c(),SampleType=c(),NF1Region=c(),Position=c(),Chromosome=c(),Patient=c())

for(i in sample.names){
    #first get general population distribution
    list.of.lists$LogRatio=c(list.of.lists$LogRatio,lrr[,i])
    list.of.lists$BAlleleFreq=c(list.of.lists$BAlleleFreq,baf[,i])
    list.of.lists$Sample=c(list.of.lists$Sample,rep(i,nrow(lrr)))
    if (i%in%c('0103','0104','0105','0106','0107','0108','0109','0110','0111','0112')){
        st<-rep('Blood',nrow(lrr))
#        sst<-rep('Blood',nrow(chr17.lrr))
    }else{
        st<-rep('Tumor',nrow(lrr))
#        sst<-rep('Tumor',nrow(chr17.lrr))
    }
    list.of.lists$SampleType=c(list.of.lists$SampleType,st)
    list.of.lists$NF1Region=c(list.of.lists$NF1Region,in.region)
    list.of.lists$Position=c(list.of.lists$Position,all.pos)
    list.of.lists$Chromosome=c(list.of.lists$Chromosome,all.chr)
    ##now get frequencies in that region of interest
 #   chr17.lists$LogRatio<-c(chr17.lists$LogRatio,chr17.lrr[,i])
 #   chr17.lists$BAlleleFreq=c(chr17.lists$BAlleleFreq,chr17.baf[,i])
 #   chr17.lists$Sample=c(chr17.lists$Sample,rep(i,nrow(chr17.lrr)))
 #   chr17.lists$SampleType=c(chr17.lists$SampleType,sst)


}
df<-data.frame(list.of.lists)
#chr17.df<-data.frame(chr17.lists)

#now plot the data
pdf('dermalNF_cnv_plots.pdf')
m<-ggplot(df,aes(x=LogRatio,colour=SampleType,linetype=NF1Region))
m<-m + geom_density() + xlim(-2.5,2)
print(m)

m<-ggplot(df,aes(x=BAlleleFreq,colour=SampleType,linetype=NF1Region))
m<-m + geom_density() + xlim(-.1,1.1)
print(m)

#now let's try histograms
m<-ggplot(df,aes(x=LogRatio,fill=SampleType,alpha=0.5,linetype=NF1Region))
m<-m + geom_histogram() + xlim(-2.5,2)
print(m)

m<-ggplot(df,aes(x=BAlleleFreq,fill=SampleType,alpha=0.5,linetype=NF1Region))
m<-m + geom_histogram() + xlim(-.1,1.1)
print(m)


dev.off()

##now focus on chr17
pdf('dermalNF_tumor_cnv_plots.pdf')
tumor.df<-subset(df,SampleType=='Tumor')

m<-ggplot(tumor.df,aes(x=LogRatio,colour=NF1Region))
m<-m + geom_density() + xlim(-2.5,2)
print(m)

m<-ggplot(tumor.df,aes(x=BAlleleFreq,colour=NF1Region))
m<-m + geom_density() + xlim(-.1,1.1)
print(m)

#now let's try histograms
m<-ggplot(tumor.df,aes(x=LogRatio,fill=NF1Region,alpha=0.5))
m<-m + geom_histogram() + xlim(-2.5,2)
print(m)

m<-ggplot(tumor.df,aes(x=BAlleleFreq,fill=NF1Region,alpha=0.5))
m<-m + geom_histogram() + xlim(-.1,1.1)
print(m)


dev.off()

pdf('dermalNF_blood_cnv_plots.pdf')
blood.df<-subset(df,SampleType=='Blood')

m<-ggplot(blood.df,aes(x=LogRatio,colour=NF1Region))
m<-m + geom_density() + xlim(-2.5,2)
print(m)

m<-ggplot(blood.df,aes(x=BAlleleFreq,colour=NF1Region))
m<-m + geom_density() + xlim(-.1,1.1)
print(m)

#now let's try histograms
m<-ggplot(blood.df,aes(x=LogRatio,fill=NF1Region,alpha=0.5))
m<-m + geom_histogram() + xlim(-2.5,2)
print(m)

m<-ggplot(blood.df,aes(x=BAlleleFreq,fill=NF1Region,alpha=0.5))
m<-m + geom_histogram() + xlim(-.1,1.1)
print(m)


dev.off()

#take most informative plot, do for each patient...

###

cna <- CNA(lrr[is.autosome,], as.character(annot$chr)[is.autosome], annot$pos[is.autosome], data.type="logratio",sample.names)
rm(annot)

smoothed.cna <- smooth.CNA(cna)
segment.smoothed.cna <- segment(smoothed.cna, verbose=1)


chr17.smoothed=subset(segment.smoothed.cna,chromlist=c("17"))
pdf('chr17.seg.smoothed.pdf')
plot(chr17.smoothed, plot.type="s")
dev.off()



segment.smoothed.cna.sundo <- segment(smoothed.cna, undo.splits="sdundo",undo.SD=2,verbose=1)
chr17.smoothed.sundo=subset(segment.smoothed.cna.sundo,chromlist=c("17"))

pdf('chr17.seg.smoothed.sundo.pdf')
plot(chr17.smoothed.sundo, plot.type="s")
dev.off()

##now get the region around chr17 to get region of interest, then re-plot b-allele frequency and logR


write.table(segment.smoothed.cna$output, file="dermal_nf_cbs_noundo.seg",
            sep="\t",quote=FALSE,row.names=F)
write.table(segment.smoothed.cna.sundo$output, file="dermal_nf_cbs_undosd2.seg",
            sep="\t",quote=FALSE,row.names=F)


## globalsd <- sapply(sample.names, function(s){
##   sd(na.omit(cna[,s]))
## })
## par(mar=c(10,4,4,4))
## barplot(globalsd,ylab="LRR standard dev",main="QC",names.arg=clnames,las=2,cex.names=.7,angle=45)

segdata <- segment.smoothed.cna$output
segdata2 <- segment.smoothed.cna.sundo$output

regions<-intersect(which(segdata2$loc.start>29000019),intersect(which(segdata2$loc.end<30427403),which(segdata2$chrom==17)))
regdata<-segdata2[regions,]

## plot(segment.smoothed.cna.sundo, plot.type="s",ylim=c(-2,2),xmaploc=TRUE)
## byChr <- subset(segment.smoothed.cna.sundo,chromlist=c("17"))
## pdf("~/foo.pdf",width=8,height=6)
## plot(byChr,plot.type="s",xmaploc=TRUE)
## dev.off()

## pdf("foo.pdf",width=10,height=8)
## for(i in 1:22){
##   byChr <- subset(segment.smoothed.cna.sundo,chromlist=c(as.character(i)))
##   plot(density(byChr$output$seg.mean),main=i)
## }
## dev.off()

#tmp <- subset(segment.smoothed.cna,samplelist=c("Sample.10"))
#plot(tmp,plot.type="w")

########################################

# cluster analysis
cnseg <- CNSeg(regdata)
rdseg <- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what = "mean")
segM <- rs(rdseg)

M <- t(do.call("rbind", lapply(4:ncol(segM), function(i) as.numeric(as.character(segM[,i])))))

idxs <- match(colnames(segM)[-1:-3], sample.names)
colnames(M) <- clnames[idxs]

plot(hclust(dist(t(M)),method="ward.D2"))

dissimilarity <- 1 - cor(M)
distance <- as.dist(dissimilarity)
plot(hclust(distance,method="ward.D2"))
