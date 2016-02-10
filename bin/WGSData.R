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
storeSomMutationFiles<-function(impact='HIGH',patient=NA){

   if(!is.na(patient))
      som.mafs=som.mafs[grep(paste('0',patient,sep=''),som.mafs$entity.patientId),]

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
#  if(exists("som.germ"))
#    return(som.germ)
  allmuts<-lapply(allsoms$entity.id,function(x) {
    print(paste('read in',x))
    read.table(synGet(x)@filePath,sep=' ',header=T,quote='"')
    })

  names(allmuts)<-allsoms$entity.id

  ##now split out somatic or germline
  som.germ<<-lapply(allmuts,function(x){
 #  print(paste('Separating out germline/somatic for sample',x))
    is.germ=apply(x,1,function(y){
      (y[['Match_Norm_Seq_Allele1']]!=y[['Reference_Allele']] || y[['Reference_Allele']]!=y[['Match_Norm_Seq_Allele2']] )})

    is.som=apply(x,1,function(y){
      (y[['Tumor_Seq_Allele1']]!=y[['Match_Norm_Seq_Allele1']] || y[['Tumor_Seq_Allele2']]!=y[['Match_Norm_Seq_Allele2']] )})
    return(list(Somatic=x[which(is.som),],Germline=x[which(is.germ),]))
  })

  return(som.germ)
}

getSomaticVars<-function(impact='HIGH',top=100){
  ##DEPRACATED
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
  df=data.frame(Gene=genes,Patient=pats,Mutation_Type=mutType,Position=pos,Tissue=tissue,MutationClass=classes)

}

getMutationsStatsByPatientAndGene<-function(patient='11',gene='NF1'){
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  print(paste('Selecting from',nrow(allsoms),'mutation files'))
  allsoms=allsoms[grep(paste('patient',patient,sep='_'),allsoms$entity.name),]
  print(paste("Found",nrow(allsoms),'for patient',patient))

  impacts<-c("LOW","MODERATE","HIGH")
  som.germ=lapply(impacts,function(x) getAllMutData(allsoms[grep(x,allsoms$entity.name),]))
  names(som.germ)=impacts

  classes=c()
  imps<-c()
  tissue=c()
  pos=c()
  ppos=c()
  mutType=c()
  dbs=c()
  exs=c()

  ##now assemble data for gene, sample
  for(i in impacts){
    for(tiss in names(som.germ[[i]])){
      ##Get somatic
      for(mt in c("Somatic","Germline")){
        gmat=som.germ[[i]][[tiss]][[mt]]
        mvals=which(gmat[,'Hugo_Symbol']==gene)
        if(length(mvals)>0){
          classes=c(classes,as.character(gmat[mvals,'Variant_Classification']))
          pos=c(pos,as.character(gmat[mvals,'HGVSc']))
          ppos=c(ppos,as.character(gmat[mvals,'HGVSp_Short']))
          exs=c(exs,as.character(gmat[mvals,'Exon_Number']))
          imps=c(imps,rep(i,length(mvals)))
          dbs<-c(dbs,as.character(gmat[mvals,'dbSNP_RS']))
          tissue=c(tissue,rep(tiss,length(mvals)))
          mutType=c(mutType,rep(mt,length(mvals)))
        }
      }
    }
  }
  df=data.frame(dbSNP=dbs,Exon=exs,Impacts=imps,Mutation_Type=mutType,Position=pos,Tissue=tissue,MutationClass=classes)

  return(df)
}


##primary function to get mutation statistics for gene
getMutationStatsForGene<-function(gene='NF1',impact=c('HIGH','MODERATE','LOW'),doPlot=TRUE,som.germ=NULL){

 if(is.null(som.germ)){
  allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
  print(paste('Selecting from',nrow(allsoms),'mutation files'))
  allsoms=allsoms[unlist(sapply(impact,grep,allsoms$entity.name)),]
  print(paste("Found",nrow(allsoms),'with',paste(impact,collapse=' or '),'impact'))
  som.germ=getAllMutData(allsoms)
  #df<-apply(allsoms,1,function(x){
 }
  classes=c()
  pats<-c()
  tissue=c()
  pos=c()
  ppos=c()
  mutType=c()
  t_depth=c()

  mut_chrom=c()
  mut_start=c()
  mut_end=c()
  ref_al=c()
  var_al=c()
  sid=c()
  for(i in 1:nrow(allsoms)){
    x=allsoms[i,]

    arr=unlist(strsplit(x[['entity.name']],split='_'))
    mv=som.germ[[x[['entity.id']]]]
    for(mt in c("Somatic","Germline")){
      mvl=mv[[mt]]
      idx=c()
      try(idx<-which(mvl[,'Hugo_Symbol']==gene))
      if(length(idx)>0){
        classes=c(classes,as.character(mvl[idx,'Variant_Classification']))
        pos=c(pos,as.character(mvl[idx,'HGVSc']))
        ppos=c(ppos,gsub('p.','',as.character(mvl[idx,'HGVSp_Short']),fixed=T))
        t_depth=c(t_depth,as.numeric(as.character(mvl[idx,'t_depth'])))
	if(mt=='Somatic')
           sid=c(sid,rep(paste(arr[1:4],collapse='_'),length(idx)))
	else
	   sid=c(sid,rep(paste(arr[1:2],collapse='_'),length(idx)))	   
        mut_chrom=c(mut_chrom,as.character(mvl[idx,'Chromosome']))
        mut_start=c(mut_start,mvl[idx,'Start_Position'])
        mut_end=c(mut_end,mvl[idx,'End_Position'])
        ra=as.character(mvl[idx,'Reference_Allele'])
        ref_al=c(ref_al,ra)
        var_al=c(var_al,apply(mvl[idx,grep('_Allele',colnames(mvl))[1:3]],1,function(y){ if(y[1]!=y[2]) return(y[2]) else return(y[3])}))

        pats=c(pats,rep(arr[2],length(idx)))
        tissue=c(tissue,rep(arr[4],length(idx)))
        mutType=c(mutType,rep(mt,length(idx)))
      }

    }
  }
  if(length(pats)==0)
    return(NULL)
  df=data.frame(Hugo_Symbol=rep(gene,length(mutType)), Protein_Change=ppos,
      Sample_ID=sid,
      Mutation_Status=mutType,Chromosome=mut_chrom,
      Start_Position=mut_start,End_Position=mut_end,
      Reference_Allele=ref_al,Variant_Allele=var_al,
      Mutation_Type=classes,TumorDepth=t_depth,
      Position=pos,Tissue=tissue,Patient=pats)

  mindf=unique(df[,-c(11,13,14)])
  write.table(mindf,file=paste(gene,paste(impact,collapse='_'),'mutations.tsv',sep=''),quote=FALSE,sep='\t',row.names=F)
  if(length(which(mutType=='Somatic'))>0){
	red.df<-subset(mindf,Mutation_Status=="Somatic")
  	write.table(red.df,file=paste(gene,paste(impact,collapse='_'),'SOMATICmutations.tsv',sep=''),quote=FALSE,sep='\t',row.names=F)
  }
  if(length(which(mutType=='Germline'))>0){
  	red.df<-subset(mindf,Mutation_Status=='Germline')
	write.table(red.df,file=paste(gene,paste(impact,collapse='_'),'GERMLINEmutations.tsv',sep=''),quote=FALSE,sep='\t',row.names=F)
  }
  df$Position=as.character(df$Position)
  mns=grep("NN+",df$Position)
  df$Position[mns]=unlist(sapply(df$Position[mns],function(x) {
    ml=attr(regexpr("N+",x),'match.length')
    gsub("N+",paste("N",ml,sep=''),x)
        }))
  ins=grep("GGTTACTCTGTTTGATTCTCGGC",df$Position)
  df$Position[ins]<-unlist(sapply(df$Position[ins],function(x) gsub("GGTTACTCTGTTTGATTCTCGGC","GGT...",x)))
  # df$Sample_ID=as.numeric(as.character(df$Sample_ID))
  require(ggplot2)
  impact=paste(impact,collapse='_or_')
  if(doPlot){
   pdf(paste('numberOf',impact,'impact',gene,'MutationsPerPatient.pdf',sep=''))
  p<-ggplot(df)+geom_bar(aes(Sample_ID,fill=Mutation_Status),position='dodge')
  p<-p+ggtitle(paste('Number of',impact,'impact  mutations in',gene))
  print(p)
  dev.off()

  ##now do class of mutation
   pdf(paste('typeOf',impact,'impact',gene,'GermlineMutationsPerPatient.pdf',sep=''))
   p<-ggplot(unique(subset(df,Mutation_Status=='Germline')))+geom_bar(aes(Sample_ID,fill=Mutation_Type),position='dodge')
   p<-p+ggtitle(paste('Type of Germline mutations in',gene))
    print(p)
    dev.off()

    pdf(paste('typeOf',impact,'impact',gene,'SomaticMutationsPerPatient.pdf',sep=''))
    p<-ggplot(subset(df,Mutation_Status=='Somatic'))+geom_bar(aes(Sample_ID,fill=Mutation_Type),position='dodge')
    p<-p+ggtitle(paste('Type of',impact,'impact Somatic mutations in',gene))
    print(p)
    dev.off()

  ##now try to classify the position of the mutation for each gene
    ##now do class of mutation

  pdf(paste('locationOf',impact,'impact',gene,'SomaticMutationsPerPatient.pdf',sep=''))
  p=ggplot(subset(df,Mutation_Status=='Somatic'))+geom_bar(aes(Sample_ID,fill=Position),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
  print(p)
  dev.off()
  pdf(paste('locationOf',impact,'impact',gene,'GermlineMutationsPerPatient.pdf',sep=''))
  p=ggplot(subset(df,Mutation_Status=='Germline'))+geom_bar(aes(Sample_ID,fill=Position),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
  print(p)
  dev.off()

  ##frequency of hits
  pdf(paste('frequencyOf',impact,'impact',gene,'SomaticMutationsPerPatient.pdf',sep=''))
  p=ggplot(subset(df,Mutation_Status=='Somatic'))+geom_bar(aes(Position,fill=Sample_ID))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
  print(p)
  dev.off()
  pdf(paste('frequencyOf',impact,'impact',gene,'GermlineMutationsPerPatient.pdf',sep=''))
  p=ggplot(subset(df,Mutation_Status=='Germline'))+geom_bar(aes(Position,fill=Sample_ID))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
  print(p)
  dev.off()

  }
  return(df)
}

heatmapFromMutDf<-function(df,fname=''){
  muts<-unique(df$Position)
  require(pheatmap)
  samp<-unique(apply(df,1,function(x) paste('Patient',paste(x[c('Patient','Tissue')],collapse='_samp_'),sep='_')))
  mmat=sapply(muts,function(x){
    mr=which(df$Position==x)
    sapply(samp,function(y){
      arr=unlist(strsplit(y,split='_'))
      rv=intersect(mr,intersect(which(df$Patient==arr[2]),which(df$Tissue==arr[4])))
    #  print(rv)
      if(length(rv)==1)
        return(as.numeric(as.character(df[rv,'TumorDepth'])))
      else
        return(0)
    })})
  ##now get the annotation for each mutation
  mut_ann=sapply(colnames(mmat),function(x) as.character(df[match(x,as.character(df$Position)),'MutationClass']))
  mmat<-mmat[order(rownames(mmat)),]
  mmat<-mmat[,order(colSums(mmat))]
  pheatmap(mmat,cluster_rows=F,cluster_cols=F,cellheight=10,cellwidth=10,annotation_col=data.frame(MutationClass=mut_ann),filename=fname)
}
