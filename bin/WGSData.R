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

getAllMutData<-function(allsoms){
  ##this function will get the number of 'high impact' somatic mutations as well
  ##as the
   #try to create matrix.
  ##first read in all files
  if(exists("som.germ"))
    return(som.germ)
  allmuts<-lapply(allsoms$entity.id,function(x) read.table(synGet(x)@filePath,sep=' ',header=T,quote='"'))
  names(allmuts)<-allsoms$entity.id
  
  ##now split out somatic or germline
  som.germ<<-lapply(allmuts,function(x){
    is.germ=apply(x,1,function(y){
      (y[['Match_Norm_Seq_Allele1']]!=y[['Reference_Allele']] || y[['Reference_Allele']]!=y[['Match_Norm_Seq_Allele2']] )})
    
    is.som=apply(x,1,function(y){
      (y[['Tumor_Seq_Allele1']]!=y[['Match_Norm_Seq_Allele1']] || y[['Tumor_Seq_Allele2']]!=y[['Match_Norm_Seq_Allele2']] )})
    return(list(Somatic=x[which(is.som),],Germline=x[which(is.germ),]))
  })
  
  return(som.germ)
}

getSomaticVars<-function(impact='HIGH',top=100){
  
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  allsoms=allsoms[grep(impact,allsoms$entity.name),]
  
  som.germ=getAllMutData(allsoms)
  classes=c()
  pats<-c()
  tissue=c()
  pos=c()
  ppos=c()
  mutType=c() 
  genes=c()
  for(i in 1:nrow(allsoms)){
    x=allsoms[i,]
    
    
    arr=unlist(strsplit(x[['entity.name']],split='_'))
    mv=som.germ[[x[['entity.id']]]]
    #first add somatic
    genes=c(genes,as.character(mv$Somatic[,'Hugo_Symbol']))
    classes=c(classes,as.character(mv$Somatic[,'Variant_Classification']))
    pos=c(pos,as.character(mv$Somatic[,'HGVSc']))
    ppos=c(ppos,as.character(mv$Somatic[,'HGVSp_Short']))
    pats=c(pats,rep(arr[2],nrow(mv$Somatic)))
    tissue=c(tissue,rep(arr[4],nrow(mv$Somatic)))
    mutType=c(mutType,rep('Somatic',nrow(mv$Somatic)))
    #then germline
    genes=c(genes,as.character(mv$Germline[,'Hugo_Symbol']))
    classes=c(classes,as.character(mv$Germline[,'Variant_Classification']))
    pos=c(pos,as.character(mv$Germline[,'HGVSc']))
    ppos=c(ppos,as.character(mv$Germline[,'HGVSp_Short']))
    pats=c(pats,rep(arr[2],nrow(mv$Germline)))
    tissue=c(tissue,rep(arr[4],nrow(mv$Germline)))
    mutType=c(mutType,rep('Germline',nrow(mv$Germline)))
    
     
  }
  df=data.frame(Gene=genes,Patient=pats,MutationType=mutType,Position=pos,Tissue=tissue,MutationClass=classes) 
  
}

getMutationStatsForGene<-function(gene='NF1',impact='HIGH'){
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  allsoms=allsoms[grep(impact,allsoms$entity.name),]
  
  som.germ=getAllMutData(allsoms)
  #df<-apply(allsoms,1,function(x){
  
  classes=c()
  pats<-c()
  tissue=c()
  pos=c()
  ppos=c()
  mutType=c()
  for(i in 1:nrow(allsoms)){
    x=allsoms[i,]
  
    arr=unlist(strsplit(x[['entity.name']],split='_'))
    mv=som.germ[[x[['entity.id']]]]
    idx.som=which(mv$Somatic[,'Hugo_Symbol']==gene)
    idx.germ=which(mv$Germline[,'Hugo_Symbol']==gene)
    print(paste('Found',length(idx.som),'somatic and',length(idx.germ),'germline mutations in patient',arr[2]))
    if(length(idx.som)>0){
      classes=c(classes,as.character(mv$Somatic[idx.som,'Variant_Classification']))
      pos=c(pos,as.character(mv$Somatic[idx.som,'HGVSc']))
      ppos=c(ppos,as.character(mv$Somatic[idx.som,'HGVSp_Short']))
      pats=c(pats,rep(arr[2],length(idx.som)))
      tissue=c(tissue,rep(arr[4],length(idx.som)))
      mutType=c(mutType,rep('Somatic',length(idx.som)))
    }
    if(length(idx.germ)>0){
      classes=c(classes,as.character(mv$Germline[idx.germ,'Variant_Classification']))
      pos=c(pos,as.character(mv$Germline[idx.germ,'HGVSc']))
      ppos=c(ppos,as.character(mv$Germline[idx.germ,'HGVSp_Short']))
      
      pats=c(pats,rep(arr[2],length(idx.germ)))
      tissue=c(tissue,rep(arr[4],length(idx.germ)))
      mutType=c(mutType,rep('Germline',length(idx.germ)))
    }
  }
  if(length(pats)==0)
    return(NULL)
  df=data.frame(Patient=pats,MutationType=mutType,Position=pos,Tissue=tissue,MutationClass=classes) 
  
  df$Position=as.character(df$Position)
  mns=grep("NN+",df$Position)
  df$Position[mns]=unlist(sapply(df$Position[mns],function(x) {
    ml=attr(regexpr("N+",x),'match.length')
    gsub("N+",paste("N",ml,sep=''),x)
        }))
  ins=grep("GGTTACTCTGTTTGATTCTCGGC",df$Position)
  df$Position[ins]<-unlist(sapply(df$Position[ins],function(x) gsub("GGTTACTCTGTTTGATTCTCGGC","GGT...",x)))
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
 
    png(paste('typeOf',gene,'SomaticMutationsPerPatient.png',sep=''))
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
  png(paste('locationOf',gene,'GermlineMutationsPerPatient.png',sep=''))
  p=ggplot(subset(df,MutationType=='Germline'))+geom_bar(aes(Patient,fill=Position),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
  print(p)
  dev.off()
  
  ##frequency of hits
  png(paste('frequencyOf',gene,'SomaticMutationsPerPatient.png',sep=''))
  p=ggplot(subset(df,MutationType=='Somatic'))+geom_bar(aes(Position,fill=Patient))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
  print(p)
  dev.off()
  png(paste('frequencyOf',gene,'GermlineMutationsPerPatient.png',sep=''))
  p=ggplot(subset(df,MutationType=='Germline'))+geom_bar(aes(Position,fill=Patient))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
  print(p)
  dev.off()
  
  return(df)
}


