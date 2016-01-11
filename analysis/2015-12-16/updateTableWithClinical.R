##take patient data and add to table
library(data.table)
library(synapseClient)

synapseLogin()


#samp.tab<-'../../data/dermalNFpatientData_addition 12-14-15.txt'
#newf=File(samp.tab,parentId='syn4984723')
#res=storeEntity(newf)

#uncomment if i delete old table before
#source("../../bin/crossDataMapping.R")
#storePatientMappingDataInTable()

clin.vars<-as.data.frame(fread(synGet('syn5555588')@filePath,strip.white=TRUE))

##now download existing data
schemaId='syn5556205'
queryResult<-synTableQuery(sprintf("select * from %s", schemaId), loadResult=TRUE)
orig.df<-queryResult@values

new.tab<-apply(clin.vars,1,function(x){
    newx=as.list(c(x,RNASeq=NA,WGS=NA,Proteomics=NA,SNPArray=NA))
    x=as.list(x)
    
    if(!is.na(x$DnaID)){
      m=subset(orig.df,as.numeric(patient)==as.numeric(x$Patient) & as.numeric(DnaID)==as.numeric(x$DnaID))
      if(nrow(m)>0){
        #print(m)
        newx$WGS=m$WGS
        newx$Proteomics=m$Prot
        newx$SNPArray=m$CNV
      }
    }
    if(!is.na(x$RnaID)){
      m=subset(orig.df,as.numeric(patient)==as.numeric(x$Patient) & as.numeric(RnaID)==as.numeric(x$RnaID))
      if(nrow(m)>0){
        print(m)        
        newx$RNASeq=m$RNASeq
      }  
    }
    newx$Patient=as.numeric(newx$Patient)
    return(unlist(newx))
})

df=as.data.frame(t(new.tab))

##now add back in PSBC
bs<-subset(orig.df,DnaID=='PBMC')
patient.ids<-1:12
tum.ids<-13:16
mf<-as.matrix(df)
for(i in 1:nrow(bs)){
  br=bs[i,]
  #get patient data
  pvals=mf[which(mf[,1]==br$patient),][1,patient.ids]
  #populate NAs for tumor data
  svals<-rep(NA,length(tum.ids))
  names(svals)<-colnames(df)[tum.ids]
  #now manually add back blood info
  rvals<-list(usedforDNA="Y",usedforRNA="not used",DnaID="PBMC",RnaID=NA,RNASeq=NA,WGS=br$WGS,Proteomics=br$Prot,SNPArray=br$CNV)
  mf<-rbind(mf,unlist(c(pvals,svals,rvals)))
}
newdf<-as.data.frame(mf)



##lastly, make sure all columns have appropriate values. 

#changed 'usedfor' columns to binary
ud=rep(FALSE,nrow(newdf))
ud[which(newdf$usedforDNA=='Y')]<-TRUE
newdf$usedforDNA=ud

ur<-rep(FALSE,nrow(newdf))
ur[which(newdf$usedforRNA=='Y')]<-TRUE
newdf$usedforRNA=ur

##also update inherited, pain, itching
v=rep(FALSE,nrow(newdf))
v[which(newdf$Inherited=='Yes')]<-TRUE
newdf$Inherited=v

v=rep(FALSE,nrow(newdf))
v[which(newdf$Pain=='Yes')]<-TRUE
newdf$Pain=v

v=rep(FALSE,nrow(newdf))
v[which(newdf$Itching=='Yes')]<-TRUE
newdf$Itching=v

colnames(newdf)[11]='NumberOfPlexiforms'
colnames(newdf)[13]='TumorNumber'
colnames(newdf)[14]='TumorLocation'
colnames(newdf)[15]='Length_in_mm'


#now table is ready to upload!!
tcresult<-as.tableColumns(newdf)
cols<-tcresult$tableColumns
fileHandleId<-tcresult$fileHandleId
projectId<-"syn4984604"
schema<-TableSchema(name="Dermal NF Sample by Patient", parent=projectId, columns =cols)
table<-Table(schema, fileHandleId)
table<-synStore(table, retrieveData=TRUE)
