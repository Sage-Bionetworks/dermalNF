##
##Compute tumor/normal CNV data
##
##
library(ASCAT)
library(data.table)
source("../../bin/dermalNFData.R")


deets=synTableQuery("SELECT Patient,DnaID,SNPArray FROM syn5556216")@values
##now get data files
all.files<-cnv_unprocessed_files()

writeAscatFiles<-function(){
  #create list of tumor files tumor logr, tumor baf, normal logr, normal baf
  plist<-NULL
  for(pat in unique(deets$Patient)){
      pdeets<-subset(deets,Patient==pat)
      norm=which(pdeets$DnaID=='PBMC')
      if(length(norm)>0){
        normsyn=pdeets$SNPArray[norm]
      }else{
        print(paste('No normal for patient',pat))
        next
      }
      normfile=all.files[[normsyn]]
      normtab<-as.data.frame(fread(normfile))
      ##normal files
      nfr=paste(normsyn,'pat',pat,'normal_logr.txt',sep='_')
      nfb=paste(normsyn,'pat',pat,'normal_baf.txt',sep='_')
      
      write.table(normtab[,c(1,3)],file=nfr,sep='\t')
      write.table(normtab[,c(1,4)],file=nfb,sep='\t')
      
      tums<-deets$SNPArray[which(!is.na(as.numeric(pdeets$DnaID)))]
      tums=tums[!is.na(tums)]
      #print(tums)
      for(tu in tums){
        tutab<-as.data.frame(fread(all.files[[tu]]))
        tur=paste(tu,'pat',pat,'tumor_logr.txt',sep='_')
        tub=paste(tu,'pat',pat,'tumor_baf.txt',sep='_')
        write.table(tutab[,c(1,3)],file=tur,sep='\t')
        write.table(tutab[,c(1,4)],file=tub,sep='\t')
        plist=rbind(plist,list(Patient=pat,TuRat=tur,TuBaf=tub,NoRat=nfr,NoBaf=nfb))
              #print(paste('patient',pat,'norm:',all.files[[normsyn]],'tum:',all.files[[tu]]))
      }
  }
  return(plist)
}




