##first get proteomics statistics

source("../../bin/dermalNFData.R")
library(WGCNA)#for clustering
require(ggplot2)


#annotes=proteomics_annotations()

#
full.df=prot_unnormalized()

#we can try clustering
#require(reshape2)
#rat.mat<-acast(full.df,Protein~Sample,value.var='Ratio')
#raw.mat<-acast(full.df,Protein~Sample,value.var='RawValue')
require(dplyr)


#zprots=which(apply(raw.mat,1,function(x) all(x==0)))

meanVals=full.df %>% group_by(Protein) %>% summarize(m=mean(RawValue))
zprots=meanVals$Protein[which(meanVals$m==0)]

min.df=full.df[-which(full.df$Protein%in%names(zprots)),]

##get variability within patient samples and across patient samples
##before and after normalization

for(metric in c("coeffOfVar","var","sd")){
  if(metric=='coeffOfVar')
    coeff<-function(x) sd(x)/mean(x)
  else if(metric=='var')
    coeff<-function(x) var(x)
  else
    coeff<-function(x) sd(x)
  
  varPatientStats=full.df %>% group_by(Patient) %>% summarize(raw=coeff(RawValue),ratio=coeff(Ratio))
  varExperimentStats=full.df %>% group_by(Experiment) %>% summarize(raw=coeff(RawValue),ratio=coeff(Ratio))

  ##now reshape into single table
  res1=tidyr::gather(varPatientStats,"ValueType","Value",2:3)
  res2=tidyr::gather(varExperimentStats,"ValueType","Value",2:3)
  colnames(res1)[1]<-colnames(res2)[1]<-'Grouping'
  full.res=rbind(cbind(Group=rep('Patient',nrow(res1)),res1),
               cbind(Group=rep("Experiment",nrow(res2)),res2))
  ratio.p=t.test(subset(res1,ValueType=='ratio')$Value,subset(res2,ValueType=='ratio')$Value)$p.value
  raw.p=t.test(subset(res1,ValueType=='raw')$Value,subset(res2,ValueType=='raw')$Value)$p.value
  pat.p=t.test(subset(res1,ValueType=='ratio')$Value,subset(res1,ValueType=='raw')$Value)$p.value
  exp.p=t.test(subset(res2,ValueType=='ratio')$Value,subset(res2,ValueType=='raw')$Value)$p.value
  
  newdf=data.frame(ValueType=c("ratio",'raw'),
                   PValue=paste("P =",format(c(ratio.p,raw.p),ndigits=4)))
  otherdf=data.frame(Group=c('Patient','Experiment'),
                     PValue=paste("P =",format(c(pat.p,exp.p),ndigits=4)))
  
  p<-ggplot(full.res)+geom_boxplot(aes(x=ValueType,fill=Group,y=Value))
  p<-p+ggtitle(paste(metric,'variability in protein measurements'))#+geom_text(aes(x=ValueType,y=PValue))
  p<-p+geom_text(data=newdf,aes(x=ValueType,y=8,label=PValue))
 
  q<-ggplot(full.res)+geom_boxplot(aes(x=Group,fill=ValueType,y=Value))
  q<-q+ggtitle(paste(metric,'variability in protein measurements'))#+geom_text(aes(x=ValueType,y=PValue))
  q<-q+geom_text(data=otherdf,aes(x=Group,y=8,label=PValue))
  
  pdf(paste(metric,'AcrossPatientsVsExpts.pdf',sep=''))
  print(p)
  print(q)
  dev.off()
}


##now update CNV analysis
###now reformulate CNV as boxplot, order by patient