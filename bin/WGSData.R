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
  
  ##now split out somatic or germline
  som.germ<-lapply(allmuts,function(x){
    is.germ=apply(x,1,function(y){
      (y[['Tumor_Seq_Allele1']]==y[['Match_Norm_Seq_Allele1']] && y[['Tumor_Seq_Allele2']]==y[['Match_Norm_Seq_Allele2']] )})
    return(list(Somatic=x[-which(is.germ),],Germline=x[which(is.germ),]))
    })
  

  
  #df<-apply(allsoms,1,function(x){
  classes=c()
  pats<-c()
  tissue=c()
  pos=c()
  mutType=c()
  for(i in 1:nrow(allsoms)){
    x=allsoms[i,]
  
    arr=unlist(strsplit(x[['entity.name']],split='_'))
    mv=som.germ[[x[['entity.id']]]]
    idx.som=which(mv$Somatic[,'Hugo_Symbol']==gene)
    idx.germ=which(mv$Germline[,'Hugo_Symbol']==gene)
    if(length(idx.som)>0){
      classes=c(classes,as.character(mv$Somatic[idx.som,'Consequence']))
      pos=c(pos,as.character(mv$Somatic[idx.som,'HGVSc']))
      pats=c(pats,rep(arr[2],length(idx.som)))
      tissue=c(tissue,rep(arr[4],length(idx.som)))
      mutType=c(mutType,rep('Somatic',length(idx.som)))
    }
    if(length(idx.germ)>0){
      classes=c(classes,as.character(mv$Germline[idx.germ,'Consequence']))
      pos=c(pos,as.character(mv$Germline[idx.germ,'HGVSc']))
      pats=c(pats,rep(arr[2],length(idx.germ)))
      tissue=c(tissue,rep(arr[4],length(idx.germ)))
      mutType=c(mutType,rep('Germline',length(idx.germ)))
    }
  }
  
  df=data.frame(Patient=pats,MutationType=mutType,Position=pos,Tissue=tissue,MutationClass=classes) 
  
  
 # df$Patient=as.numeric(as.character(df$Patient))
  require(ggplot2)
  png(paste('numberOf',gene,'MutationsPerPatient.png',sep=''))
  p<-ggplot(df)+geom_bar(aes(Patient,fill=MutationType),position='dodge')
  p<-p+ggtitle(paste('Number of high-impact  mutations in',gene))
  print(p)
  dev.off()

  ##now do class of mutation
   png(paste('typeOf',gene,'GermlineMutationsPerPatient.png',sep=''))
   p<-ggplot(unique(subset(df,MutationType=='Germline')))+geom_bar(aes(Patient,fill=MutationClass),position='dodge')
   p<-p+ggtitle(paste('Type of high-impact Germline mutations in',gene))
    print(p)
    dev.off()
 
    png(paste('typeOf',gene,'SomaticMutationsPerPatient.png',sep=''),width=1000)
    p<-ggplot(subset(df,MutationType=='Somatic'))+geom_bar(aes(Patient,fill=MutationClass),position='dodge')
    p<-p+ggtitle(paste('Type of high-impact Somatic mutations in',gene))
    print(p)
    dev.off()   
    
  ##now try to classify the position of the mutation for each gene
    ##now do class of mutation

  png(paste('locationOf',gene,'SomaticMutationsPerPatient.png',sep=''))
  p=ggplot(subset(df,MutationType=='Somatic'))+geom_bar(aes(Patient,fill=Position),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
  print(p)
  dev.off()
  png(paste('locationOf',gene,'GermlineMutationsPerPatient.png',sep=''),width=1000)
  p=ggplot(subset(df,MutationType=='Germline'))+geom_bar(aes(Patient,fill=Position),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
  print(p)
  dev.off()
  return(df)
}


