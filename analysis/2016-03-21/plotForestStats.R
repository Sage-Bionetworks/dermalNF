##evaluate forest parameters

library(synapseClient)
synapseLogin()
library(dplyr)
library(ggplot2)

stats<- read.table(synGet('syn5816309')@filePath,sep='\t',header=T,fill=T)
stats$params=apply(stats,1,function(x) paste(x[3:5],collapse='_'))

  ##we can't really compare networks across samples to across patients
acrossSamples=stats[grep('AcrossSamples',stats$patientCoverage),]
perPatient<-stats[grep('AcrossPatient',stats$patientCoverage),]

pdf('edgeAndTreeSizesByparam.pdf')
p<-ggplot(acrossSamples)+geom_point(aes(y=numEdges,x=params,color=inputDataType,shape=hasUBC))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p+ggtitle('Num Edges Across Samples'))
p<-ggplot(acrossSamples)+geom_point(aes(x=params,y=numTrees,color=inputDataType,shape=hasUBC))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p+ggtitle('Num Trees Across Samples'))
p<-ggplot(perPatient)+geom_point(aes(y=numEdges,x=params,color=inputDataType,shape=hasUBC))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p+ggtitle('Num Edges per patient'))
p<-ggplot(perPatient)+geom_jitter(aes(x=params,y=numTrees,color=inputDataType,shape=hasUBC),height=0.1,width=0.3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p+ggtitle('num Trees per patient'))
dev.off()

##get original files
in.files=list.files("../2016-03-18")
in.files=in.files[grep('Proteins',in.files)]
mets=c()
vals=c()
groups=c()
for(a in in.files[grep('.tab',in.files)]){
  val=read.table(paste('../2016-03-18',a,sep='/'),sep='\t')[,2]
  met=unlist(strsplit(basename(a),split='OfProteins'))[1]
  if(length(grep("Patient",a))>0){
    groups=c(groups,rep('Patient',length(val)))
  }else{
    groups=c(groups,rep('All',length(val)))
  }
  mets=c(mets,rep(met,length(val)))
  vals<-c(vals,val)
}
input.data<-data.frame(Group=groups,Metric=mets,Values=vals)

pdf('distributionOfOriginalValues.pdf')
p<- ggplot(subset(input.data,Metric=='sum'))+geom_histogram(aes(Values,fill=Group))+ggtitle('Distribution of sum scores')
print(p)
p<- ggplot(subset(input.data,Metric=='mean'))+geom_histogram(aes(Values,fill=Group))+ggtitle('Distribution of mean scores')
print(p)
p<- ggplot(subset(input.data,Metric=='frac'))+geom_histogram(aes(Values,fill=Group))+ggtitle('Distribution of fractional scores')
print(p)

dev.off()


