source("../../bin/crossDatasetAnalysis.R")

library(ggplot2)
##compare and plot distribtuions with significnace
compareDists<-function(orig,rand,datatype=''){
  #first do patient
  patdf<-data.frame(Sample=c(rep('Original',length(orig$patientCorrelation)),
                             rep('Random',length(rand$patientCors))),
                    Correlation=c(orig$patientCorrelation,rand$patientCors))
  pval=t.test(orig$patientCorrelation,rand$patientCors)$p.value
  
  png(paste(datatype,'perPatientCors.png',sep=''))
  p<-ggplot(patdf) +
    geom_histogram(aes(x=Correlation,fill=Sample),alpha=0.5,stat='density') +  
    ggtitle(paste('Per patient data correlation compared to random,\np =',format(pval,ndigits=3)))
  print(p)
  dev.off()
  
  #now do gene
  genedf<-data.frame(Sample=c(rep('Original',length(orig$geneCorrelation)),
                             rep('Random',length(rand$geneCors))),
                    Correlation=c(orig$geneCorrelation,rand$geneCors))
  pval=t.test(orig$geneCorrelation,rand$geneCors)$p.value
  
  png(paste(datatype,'perGeneCors.png',sep=''))
  p<-ggplot(genedf) +
    geom_histogram(aes(x=Correlation,fill=Sample),alpha=0.5,stat='density') +  
    ggtitle(paste('Per gene data correlation compared to random,\np =',format(pval,ndigits=3)))
  print(p)
  dev.off()
}

##1 execute permutation test to determine if correlation distribution is significant...

dat.mat<-matchMrnaProt()
mat1<-dat.mat$rna
mat2<-dat.mat$prot

rand.dists<-getMatchedDistributions(mat1,mat2,numiter=1000) ##this will return gene and patient correlations
#2 which genes are more correlated than we'd expect by chance?  
#now get original correlations
mp.orig.dists<-mrna_prot_comp(doPlot=FALSE)

compareDists(mp.orig.dists,rand.dists,'mRNAProtein')


##now move to gene and CNV
dat.mat<-matchMrnaCnv()
mat1<-dat.mat$rna
mat2<-dat.mat$cnv

rand.dists<-getMatchedDistributions(mat1,mat2,numiter=1000) ##this will return gene and patient correlations
#2 which genes are more correlated than we'd expect by chance?  
#now get original correlations
mp.orig.dists<-mrna_prot_comp(doPlot=FALSE)

compareDists(mp.orig.dists,rand.dists,'mRNACNV')

#3 which genes are less correlated? 

##4 linear model of RNA-Seq variation