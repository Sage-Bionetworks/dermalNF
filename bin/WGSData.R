##process WGS data

library(synapseClient)
synapseLogin()

allmafs<-synapseQuery("select * from entity where parentId=='syn5522808'")
som.mafs<-allmafs[which(allmafs$entity.tissueType=='tumorVsNormal'),]
gl.mafs<-allmafs[which(allmafs$entity.tissueType=='PBMC'),]


getMutationSummary<-function(){

    allMuts<-sapply(allmafs$entity.id,function(x){
        res<-synGet(x)
        tab<-read.table(gzfile(res@filePath),sep='\t',header=T)

    })
}
#summary(tab$Consequence)

library(parallel)
library(data.table)
storeSomMutationFiles<-function(impact='HIGH'){

    allMuts<-lapply(som.mafs$entity.id,function(x){
        res<-synGet(x)
       # if(res@annotations$tissueType=='PBMC')
      #      return(NULL)
        fp=res@filePath
        fname=paste('patient',gsub('CT0+','',res@annotations$patientId),'tissue',res@annotations$tissueID,impact,'impact_somaticMutations.maf',sep='_')
        if(! file.exists(fname)){
            tab<-as.data.frame(fread(paste('zcat',fp)))
                                        #dont filter by consequence
                                        #vars<-tab[which(tab$Consequence%in%mutClasses),]
            if(!is.na(impact))
                vars<-tab[which(tab$IMPACT==impact),]
            else
                vars<-tab
            print(length(which(vars$Hugo_Symbol=='NF1')))#summary(vars$Hugo_Symbol))

            write.table(vars,file=fname,row.names=F,quote=F)
        }
        sf=File(fname,parentId='syn5578958')
        annotations(sf)<-res@annotations
        executed(sf)<-'https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/WGSData.R'
        synStore(sf)
        return (vars)
    })#,mc.cores=4)
    return(allMuts)
}

##now do the same for the germline files
storeGermlineMutationFiles<-function(impact='HIGH'){
  
  allMuts<-lapply(gl.mafs$entity.id,function(x){
    res<-synGet(x)
    # if(res@annotations$tissueType=='PBMC')
    #      return(NULL)
    fp=res@filePath
    fname=paste('patient',gsub('CT0+','',res@annotations$patientId),'tissue',res@annotations$tissueID,impact,'impact_germlineMutations.maf',sep='_')
    if(! file.exists(fname)){
      tab<-as.data.frame(fread(paste('zcat',fp)))
      #dont filter by consequence
      #vars<-tab[which(tab$Consequence%in%mutClasses),]
      if(!is.na(impact))
        vars<-tab[which(tab$IMPACT==impact),]
      else
        vars<-tab
      print(length(which(vars$Hugo_Symbol=='NF1')))#summary(vars$Hugo_Symbol))
      
      write.table(vars,file=fname,row.names=F,quote=F)
    }
    sf=File(fname,parentId='syn5580983')
    annotations(sf)<-res@annotations
    executed(sf)<-'https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/WGSData.R'
    synStore(sf)
    return (vars)
  })#,mc.cores=4)
  return(allMuts)
}

getMutationStatsForGene<-function(gene='NF1'){
  ##this function will get the number of 'high impact' somatic mutations as well
  ##as the
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  #try to create matrix.
  ##first read in all files
  allmuts<-lapply(allsoms$entity.id,function(x) read.table(synGet(x)@filePath,sep=' ',header=T,quote='"'))
  names(allmuts)<-allsoms$entity.id

  df<-apply(allsoms,1,function(x){
    arr=unlist(strsplit(x[['entity.name']],split='_'))
    mcounts=length(which(allmuts[[x[['entity.id']]]][,'Hugo_Symbol']==gene))
    c(NumMutations=mcounts,Gene=gene,Patient=arr[2],Tissue=arr[4])
  })
  df<-as.data.frame(t(df))
  df$NumMutations=as.numeric(as.character(df$NumMutations))
 # df$Patient=as.numeric(as.character(df$Patient))
  require(ggplot2)
  png(paste('numberOf',gene,'MutationsPerPatient.png',sep=''))
  p<-ggplot(df)+geom_jitter(aes(Patient,NumMutations,colour=Patient),height=0.1)
  p<-p+ggtitle(paste('Number of high-impact somatic mutations in',gene))
  print(p)
  dev.off()

  ##now do class of mutation
  classes=c()
  pats<-c()
  tissue=c()
  for(i in 1:nrow(allsoms)){
    x=allsoms[i,]
    arr=unlist(strsplit(x[['entity.name']],split='_'))
    idx=which(allmuts[[x[['entity.id']]]][,'Hugo_Symbol']==gene)
    if(length(idx)>0){
      classes=c(classes,as.character(allmuts[[x[['entity.id']]]][idx,'Consequence']))
      pats=c(pats,rep(arr[2],length(idx)))
      tissue=c(tissue,rep(arr[4],length(idx)))
    }
  }
   newdf<-data.frame(VarType=classes,Patient=pats,Tissue=tissue)
   png(paste('typeOf',gene,'MutationsPerPatient.png',sep=''))

   p<-ggplot(newdf)+geom_bar(aes(Patient,fill=VarType),position='dodge')
   p<-p+ggtitle(paste('Type of high-impact somatic mutations in',gene))
    print(p)
    dev.off()
  return(df)
}
