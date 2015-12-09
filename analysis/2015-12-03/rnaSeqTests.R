source("../../bin/crossDatasetAnalysis.R")

library(ggplot2)
##compare and plot distribtuions with significnace
compareDists<-function(orig,rand,datatype=''){
  #first do patient
  patdf<-data.frame(Sample=c(rep('Original',length(orig$patientCorrelation)),
                             rep('Random',length(rand$patientCors))),
                    Correlation=c(orig$patientCorrelation,rand$patientCors))
  pval=wilcox.test(orig$patientCorrelation,rand$patientCors,alt='g')$p.value

  png(paste(datatype,'perPatientCors.png',sep=''))
  p<-ggplot(patdf) +
    geom_histogram(aes(x=Correlation,fill=Sample),stat='density') +
    ggtitle(paste('Per sample data correlation compared to random,\np =',format(pval,ndigits=3)))
  print(p)
  dev.off()

  #now do gene
  genedf<-data.frame(Sample=c(rep('Original',length(orig$geneCorrelation)),
                             rep('Random',length(rand$geneCors))),
                    Correlation=c(orig$geneCorrelation,rand$geneCors))
  pval=wilcox.test(orig$geneCorrelation,rand$geneCors,alt='g')$p.value

  png(paste(datatype,'perGeneCors.png',sep=''))
  p<-ggplot(genedf) +
    geom_histogram(aes(x=Correlation,fill=Sample),stat='density') +
    ggtitle(paste('Per gene data correlation compared to random,\np =',format(pval,ndigits=3)))
  print(p)
  dev.off()
}

##1 execute permutation test to determine if correlation distribution is significant...
##first for protein/mrna
reps=1000
dat.mat<-matchMrnaProt()
mat1<-dat.mat$rna
mat2<-dat.mat$prot

rand.dists<-getMatchedDistributions(mat1,mat2,numiter=reps) ##this will return gene and patient correlations
#2 which genes are more correlated than we'd expect by chance?
#now get original correlations
mp.orig.dists<-mrna_prot_comp(doPlot=FALSE)

compareDists(mp.orig.dists,rand.dists,'mRNAProtein')

###2 now find thresholds
geneThresh=sort(rand.dists$geneCors,dec=T)[0.05*reps]
patThresh=sort(rand.dists$patientCors,dec=T)[0.05*reps]

sig.genes<-names(mp.orig.dists$geneCorrelation)[which(mp.orig.dists$geneCorrelation>geneThresh)]
sig.pats<-names(mp.orig.dists$patientCorrelation)[which(mp.orig.dists$patientCorrelation>patThresh)]

pdf('genesAndProtsMoreCorrelatedByChance.pdf')
for(g in sig.genes){
    plot(mat1[,g],mat2[,g],xlab=paste(g,'mRNA Expression'),ylab=paste(g,'Protein Expression'))
}
dev.off()

##now move to gene and CNV
#first get matched data
dat.mat<-matchMrnaCnv()
mat1<-dat.mat$rna
mat2<-dat.mat$cnv

rand.dists<-getMatchedDistributions(mat1,mat2,numiter=reps) ##this will return gene and patient correlations
#2 which genes are more correlated than we'd expect by chance?
#now get original correlations
mc.orig.dists<-mrna_cnv_comp(doPlot=FALSE)

compareDists(mc.orig.dists,rand.dists,'mRNACNV')

geneThresh=sort(rand.dists$geneCors,dec=T)[0.01*reps]
patThresh=sort(rand.dists$patientCors,dec=T)[0.05*reps]

sig.genes<-names(mc.orig.dists$geneCorrelation)[which(mc.orig.dists$geneCorrelation>geneThresh)]
sig.pats<-names(mc.orig.dists$patientCorrelation)[which(mc.orig.dists$patientCorrelation>patThresh)]

pdf('genesAndCNVMoreCorrelatedByChance.pdf')
for(g in sig.genes[order(mc.orig.dists$geneCorrelation[sig.genes],decreasing=T)]){
  plot(mat1[,g],mat2[,g],xlab=paste(g,'mRNA Expression'),ylab=paste(g,'Copy number log R ratio'))
}
dev.off()

negThresh=sort(rand.dists$geneCors,dec=F)[0.01*reps]
neg.genes<-names(mc.orig.dists$geneCorrelation)[which(mc.orig.dists$geneCorrelation<negThresh)]
