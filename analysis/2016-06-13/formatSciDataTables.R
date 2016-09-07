##format tables for Sci Data paper

library(synapseClient)
synapseLogin()

dna.stats<-read.csv(synGet('syn6156139')@filePath,as.is=T)
rna.stats<-read.csv(synGet('syn6156141')@filePath)
rna.mets<-read.csv(synGet('syn6156140')@filePath)
proper.dna.ids<-read.table(synGet('syn4999547')@filePath,as.is=T,sep='\t',header=T)

dna.stats$UpdatedPatid<-proper.dna.ids$Patient.ID[match(dna.stats$Library.Name,proper.dna.ids$Sample.ID)]

##and the discordance analysis
disc<-read.table(synGet('syn5669991')@filePath)

##now get table data to combine all patient samples
tab<-synTableQuery('SELECT Patient,DnaID,RnaID,RNASeq,TumorNumber FROM syn5556216')@values
tab$DnaID<-sapply(tab$DnaID,trimws)

final.res<-apply(dna.stats,1,function(x) {
  pids<-unlist(strsplit(x[6],split=' '))
  if(length(pids)==1)
    pids<-unlist(strsplit(x[6],split='_'))
  pat=gsub("CT0+","",pids[1])
  samp=gsub('^0+','',pids[2])
  tsamp=tab$TumorNumber[intersect(which(tab$Patient==pat),which(tab$DnaID==samp))]
  if(samp=='PBMC')
    tsamp='PBMC'
  return(list(Patient=pat,Sample=tsamp,TotalReads=x[4],PercentAligned=x[5]))
})

dna.tab<-data.frame(do.call('rbind',final.res))
dna.tab<-dna.tab[order(as.numeric(dna.tab$Patient)),]
write.table(dna.tab,file='dnaStatsTable.csv',sep=',',quote=F,row.names=F)


rna.idx=match(rna.stats$Sample.Name,rna.mets$Library.Name)
rna.tab<-cbind(rna.stats,rna.mets[rna.idx,])

rna.files<-synQuery("select name,id from entity where parentId=='syn5493036'")
rna.files$library=sapply(rna.files$entity.name,function(x) unlist(strsplit(x,split='_'))[3])

rna.syn<-match(rna.tab$Library.Name,rna.files$library)

rna.full.tab<-cbind(rna.tab,rna.files[rna.syn,])

rna.pat<-tab[which(!is.na(tab$RNASeq)),]

rna.idx<-match(rna.pat$RNASeq,rna.full.tab$entity.id)

comp<-cbind(rna.pat,rna.full.tab[rna.idx,])

df<-data.frame(Patient=comp$Patient,TumorNumber=comp$TumorNumber,TotalRNABases=comp$PF_BASES,Aligned=comp$PF_ALIGNED_BASES,RIN=comp$RIN.Score)
tab$TumorNumber[which(tab$DnaID=='PBMC')]<-'PBMC'

#now merge everything into a single table!!
fulldf<-apply(tab,1,function(x){
  patient=trimws(x[1])
  tumor=trimws(x[5])
  print(paste(patient,tumor))
  dna<-intersect(which(dna.tab$Patient==patient),which(dna.tab$Sample==tumor))
  rna<-intersect(which(df$Patient==patient),which(df$TumorNumber==tumor))
  if(length(dna)==0){
    dna<-c(DnaReads=NA,PercentAligned=NA)
  }else{
    dna<-c(DnaReads=dna.tab$TotalReads[dna],PercentAligned=dna.tab$PercentAligned[dna])
  }
  if(length(rna)==0){
    rna<-c(RnaReads=NA,PercentRnaAligned=NA,RIN=NA)
  }else{
    rna<-c(RnaReads=df$TotalRNABases[rna],PercentRnaAligned=df$FracRnaAligned[rna],RIN=df$RIN[rna])
  }
  return(c(Patient=patient,Tumor=tumor,dna,rna))
    
})

newdf<-data.frame(do.call('rbind',fulldf))

names(newdf)[1:2]<-c("Patient","Tumor")
newdf$PercentRnaAligned<-sapply(as.numeric(newdf$PercentRnaAligned)*100,format,digits=4)


newdf$DnaReads=sapply(newdf$DnaReads,function(x) format(as.numeric((gsub(',','',x))),digits=4, scientific=T))

newdf$RnaReads=sapply(newdf$RnaReads,function(x) format(as.numeric((gsub(',','',x))),digits=4, scientific=T))

newdf<-apply(newdf,2,unlist)
write.table(newdf,file='combinedDnaRnaQualityTable.csv',sep=',',quote=T,row.names=F)


