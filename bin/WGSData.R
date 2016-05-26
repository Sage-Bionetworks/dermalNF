##General set of commands to analyze and process WGS data
##

library(synapseClient)
synapseLogin()
library(data.table)

require(dplyr)
require(reshape2)
require(pheatmap)
##read in all cancer gene mutations

#cancer.gene.muts<-read.table(synGet('syn5611520')@filePath,header=T,sep='\t')

if(!exists('all.gene.muts'))
  all.gene.muts<-read.table(synGet('syn5839666')@filePath,header=T,sep='\t')

if(!exists('cancer.gene.muts')){
 cancer.genes=unique(read.table('../../data/Census_allTue Jan 19 18-58-56 2016.csv',sep=',',header=T,quote='"')$Gene.Symbol)
 cancer.gene.muts<-subset(all.gene.muts,Hugo_Symbol%in%cancer.genes)
}
#  cancer.gene.muts<-read.table(synGet('syn5611520')@filePath,header=T,sep='\t')

#all.gene.muts<-read.table(synGet('syn5713423')@filePath,header=T,sep='\t')
doPatientHeatmap<-function(mut.tab,title,fname,minSamples=1){
  
  ##format into matrix
  mut.counts=mut.tab%>% 
    group_by(Hugo_Symbol,Sample_ID) %>% 
    summarize(Count=n()) %>%
    acast(Hugo_Symbol~Sample_ID,value.var='Count')
  
  #get extra variables, like distinct variants
  num.variants=mut.tab%>%group_by(Hugo_Symbol)%>%summarize(Variants=n_distinct(Protein_Change))
  variants=num.variants$Variants
  names(variants)=num.variants$Hugo_Symbol
  
  ##now filter by minCount
  mut.counts=mut.counts[which(apply(mut.counts,1,function(x) length(which(x>0)))>minSamples),]
  
  
  
  #get row and column order
  r.o=order(apply(mut.counts,1,function(x) length(which(x>0))))
  c.o=order(apply(mut.counts,2,function(x) length(which(x>0))))
  
  mut.counts=mut.counts[r.o,c.o]
  
  pheatmap(log10(1+mut.counts),cellwidth=10,cellheight=10,cluster_rows=F,cluster_cols=F,
           main=title,filename=fname,annotation_row=data.frame(Variants=variants))
             
}

#'Define function to plot gene mutations across patients in a heatmap
#'@param mutTable - table of mutations in cBio format
#'@param minSampls - minimum number of samples required to include a gene
#'@param notIncluded - what type of mutation to exclude?
#'
panPatientPlots<-function(mutTable=all.gene.muts,minSamples=2,notIncluded=c(),prefix=''){
  ##first remove genes
  if(length(notIncluded)>0){
    print(paste("Removing",length(which(mutTable$Mutation_Type%in%notIncluded)),'mutations that are',paste(notIncluded,collapse='or')))
    mutTable=mutTable[-which(mutTable$Mutation_Type%in%notIncluded),]  
  }
  
  #get somatic mutants, then plot
  som.muts=subset(mutTable,Mutation_Status=='Somatic')
  #now filter by minSamples value
  som.sample.count=som.muts%>%group_by(Hugo_Symbol)%>%summarize(Samples=n_distinct(Sample_ID))
  genes.to.plot=as.character(som.sample.count$Hugo_Symbol[which(som.sample.count$Samples>minSamples)])
  som.muts=som.muts[which(as.character(som.muts$Hugo_Symbol)%in%genes.to.plot),]
  
  title=paste('Number of somatic mutations in\n genes',
              ifelse(length(notIncluded)>0,paste('(not',paste(notIncluded,collapse=','),')'),''),
              'that occur in at least',minSamples,'samples')
  fname=paste(prefix,'somaticMuts_not',paste(notIncluded,collapse='_'),'minSamples',minSamples,sep='_')
  doPatientHeatmap(som.muts,title,paste(fname,'png',sep='.'))
  
  #then get germline
  g.muts=subset(mutTable,Mutation_Status=='Germline')
  title=paste('Number of somatic mutations in\n genes',
              ifelse(length(notIncluded)>0,paste('(not',paste(notIncluded,collapse=','),')'),''),
              'that occur in at least',minSamples,'samples')
  fname=paste(prefix,'somaticMuts_not',paste(notIncluded,collapse='_'),'minSamples',minSamples,sep='_')
 
  g.sample.count=g.muts%>%group_by(Hugo_Symbol)%>%summarize(Samples=n_distinct(Sample_ID))
  genes.to.plot=as.character(g.sample.count$Hugo_Symbol[which(g.sample.count$Samples>minSamples)])
  g.muts=g.muts[which(as.character(g.muts$Hugo_Symbol)%in%genes.to.plot),]
  
  title=paste('Number of germline mutations in\n genes',
              ifelse(length(notIncluded)>0,paste('(not',paste(notIncluded,collapse=','),')'),''),
              'that occur in at least',minSamples,'samples')
  fname=paste(prefix,'germlineMuts_not',paste(notIncluded,collapse='_'),'minSamples',minSamples,sep='_')
  doPatientHeatmap(g.muts,title,paste(fname,'png',sep='.'))
  
}


#' Get all MAF files that were calculated using VCF2MAF
#' @param mut.type can be 'all' (default),'somatic' or 'germline'
#' @return table describing entities of  all MAF files on synapse
getMAFs<-function(mut.type='all'){

  allmafs<-synapseQuery("select * from entity where parentId=='syn5522808'")
  som.mafs<-allmafs[which(allmafs$entity.tissueType=='tumorVsNormal'),]
  gl.mafs<-allmafs[which(allmafs$entity.tissueType=='PBMC'),]

  if(tolower(mut.type)=='all')
    return(allmafs)
  else if(tolower(mut.type)=='somatic')
    return(som.mafs)
  else if (tolower(mut.type)=='germline')
    return(gl.mafs)
  else
    print("mut.type must be either 'all','somatic',or 'germline'")
  return(NULL)
}

#'getMutationSummary opens all maf files and loads into list
#'@param allmafs - list of MAFs to summarize
#'@return list of tables
getMutationSummary<-function(allmafs=getMAFs('all')){

    allMuts<-sapply(allmafs$entity.id,function(x){
        res<-synGet(x)
        tab<-read.table(gzfile(res@filePath),sep='\t',header=T)

    })
}
#summary(tab$Consequence)

#library(parallel)


#'Takes original MAF files and filters by predicted 'impact' and/or patient.  Will include both
#'Somatic and Germline mutations, because they are both in the original MAF files for each tumor
#'@param som.mafs  MAFs to be analyzed - the somatic files include both the somatic and germline
#'@param impact can be "HIGH",'MODERATE' or 'LOW'
#'@param patient
storeSomMutationFiles<-function(som.mafs=getMAFs('somatic'),impact='HIGH',patient=NA){

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

#'getAllMutData is a function that will take a mutation maf table and filter
#'by whether or not a particular mutation is somatic or germline
#'@param allsoms is a table of all mutations of a particular gene, impact or patient
#'@return list of two tables 'Somatic' is table of somatic mutations while 'Germline' is table of germline
getAllMutData<-function(allsoms=getMAFs('all'),filter=c()){
  ##this function will get the number of 'high impact' somatic mutations as well
  ##as the
   #try to create matrix.
  ##first read in all files

  allmuts<-lapply(allsoms$entity.id,function(x) {
    print(paste('read in',x))
    read.table(synGet(x)@filePath,sep=' ',header=T,quote='"')
    })

  names(allmuts)<-allsoms$entity.id

  ##now split out somatic or germline
  som.germ<<-lapply(allmuts,function(x){
 #  print(paste('Separating out germline/somatic for sample',x))
    fout=which(as.character(x$FILTER)%in%filter)
    if(length(fout)>0){
        print(paste('Keeping',length(fout),'out of',nrow(x),'because they are',paste(filter,collapse=',')))
        x=x[fout,]
    }
    is.germ=apply(x,1,function(y){
      (y[['Match_Norm_Seq_Allele1']]!=y[['Reference_Allele']] || y[['Reference_Allele']]!=y[['Match_Norm_Seq_Allele2']] )})

    is.som=apply(x,1,function(y){
      (y[['Tumor_Seq_Allele1']]!=y[['Match_Norm_Seq_Allele1']] || y[['Tumor_Seq_Allele2']]!=y[['Match_Norm_Seq_Allele2']] )})
    return(list(Somatic=x[which(is.som),],Germline=x[which(is.germ),]))
  })

  return(som.germ)
}


#'getMutationStatsForGene obtains all mutations for a particular gene of interest across all patients
#'@param gene is the gene symbol in question
#'@param impact is a list of which mutations to include, defaults to all ('HIGH','MODERATE' and 'LOW')
#'@param doPlot: if set to true, will plot some basic statistics about where and when this mutation occurs
#'@param som.germ - the MAF file tables separated by whether or not the mutation is somatic or germline
getMutationStatsForGene<-function(gene='NF1',impact=c('HIGH','MODERATE','LOW'),doPlot=FALSE,filter=c(),som.germ=getAllMutData(filter=filter),redo=FALSE){

  ##first check to see if we have the file already on synapse
  if(gene%in%all.gene.muts$Hugo_Symbol && !redo){
    print(paste('Found gene',gene,' already processed, will analyze mutations of all impact (impact argument ignored'))
    df=subset(all.gene.muts,Hugo_Symbol==gene)
  }else{
    print(paste('No evidence of',gene,'in mutation data'))
    if(!redo){
      print('Set redo=TRUE to double check')
      return(data.frame())
    }
    if(is.null(som.germ)){
      allsoms<-synapseQuery("select * from entity where parentId=='syn5578958'")
      print(paste('Selecting from',nrow(allsoms),'mutation files'))
      allsoms=allsoms[unlist(sapply(impact,grep,allsoms$entity.name)),]
      print(paste("Found",nrow(allsoms),'with',paste(impact,collapse=' or '),'impact'))
      som.germ=getAllMutData(allsoms,filter=filter)
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
  
    
    ##the mindf files are visible via cbioportal.
    mindf=unique(df[,-c(11,13,14)])
    write.table(mindf,file=paste(gene,paste(impact,collapse='_'),'mutations.tsv',sep=''),quote=FALSE,sep='\t',row.names=F)
    #somatic only
    if(length(which(mutType=='Somatic'))>0){
  	  red.df<-subset(mindf,Mutation_Status=="Somatic")
    	write.table(red.df,file=paste(gene,paste(impact,collapse='_'),'SOMATICmutations.tsv',sep=''),quote=FALSE,sep='\t',row.names=F)
    }
    #germline only
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
    impact=paste(impact,collapse='_or_')

  }
  if(doPlot){
    require(ggplot2)
    
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
    p=ggplot(subset(df,Mutation_Status=='Somatic'))+geom_bar(aes(Sample_ID,fill=Protein_Change),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
    print(p)
    dev.off()
    pdf(paste('locationOf',impact,'impact',gene,'GermlineMutationsPerPatient.pdf',sep=''))
    p=ggplot(subset(df,Mutation_Status=='Germline'))+geom_bar(aes(Sample_ID,fill=Protein_Change),position='dodge')+ggtitle(paste(gene,'Mutation Position'))
    print(p)
    dev.off()

  ##frequency of hits
    pdf(paste('frequencyOf',impact,'impact',gene,'SomaticMutationsPerPatient.pdf',sep=''))
    p=ggplot(subset(df,Mutation_Status=='Somatic'))+geom_bar(aes(Protein_Change,fill=Sample_ID))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
    print(p)
    dev.off()
    pdf(paste('frequencyOf',impact,'impact',gene,'GermlineMutationsPerPatient.pdf',sep=''))
    p=ggplot(subset(df,Mutation_Status=='Germline'))+geom_bar(aes(Protein_Change,fill=Sample_ID))+ggtitle(paste(gene,'Mutation Position'))+theme(axis.text.x=element_text(angle = -90, hjust = 0))
    print(p)
    dev.off()
  }
  
  return(df)
}

#'Create a heatmap with PHRED scores (need to add in separately)
heatmapWithPhredScore<-function(df,fname,phredscore,cutoff=10,samp.vars=NA){
  require(reshape2)
  require(dplyr)
  require(pheatmap)
  
  ##first map df to phred
  mut.idx<-match(df$Start_Position,phredscore$Pos)
  na.vals<-which(is.na(mut.idx))
  mut.idx[na.vals]<-sapply(na.vals,function(x){
    match(df[x,'Start_Position']-1,phredscore$Pos)
  })
 df$PHRED<-phredscore$PHRED[mut.idx]
# print(head(df))
  #df$Consequence<-nf1.deets$Consequence[mut.idx]
  
  ##first start with matrix of phred scores
   ##then separate out by tumor type
  changes<-as.character(df$Protein_Change)
  types<-df$Mutation_Type[match(unique(changes),as.character(df$Protein_Change))]
  ttypes<-data.frame(Consequence=types)
  rownames(ttypes)<-unique(changes)
 # print(head(ttypes))
#  print(head(types))
  
  pmat = acast(df,Sample_ID~Protein_Change,value.var='PHRED',fill=0)
  col.ords=order(apply(pmat,2,function(x) mean(x[which(x>0)])))
  pmat<-pmat[order(rowSums(pmat)),col.ords]
  col.cuts=which(apply(pmat,2,function(x) mean(x[which(x>0)])>cutoff))
  
  soms=grep('tissue',rownames(pmat))
  #print(head(pmat))
  
  pheatmap(pmat[soms,],annotation_col = ttypes,annotation_row=samp.vars,
           cellheight=10,cellwidth=10,cluster_rows=F,cluster_cols=F,
         #  legend_labels=as.character(types),legend_breaks=c(type.nums)-0.5,
           filename=paste('positionScoreSomatic',fname,sep=''))
  
  pheatmap(pmat[-soms,],annotation_col = ttypes,annotation_row=samp.vars,
           cluster_rows=F,cluster_cols=F,
           cellheight=10,cellwidth=10,
         filename=paste('positionScoreGermline',fname,sep=''))
  
  pheatmap(pmat[soms,col.cuts],annotation_col = ttypes,annotation_row=samp.vars,
           cellheight=10,cellwidth=10,cluster_rows=F,cluster_cols=F,
           #  legend_labels=as.character(types),legend_breaks=c(type.nums)-0.5,
           filename=paste('positionScoreOver',cutoff,'Somatic',fname,sep=''))
  
  pheatmap(pmat[-soms,col.cuts],annotation_col = ttypes,annotation_row=samp.vars,
           cluster_rows=F,cluster_cols=F,
           cellheight=10,cellwidth=10,
           filename=paste('positionScoreOver',cutoff,'Germline',fname,sep=''))
  
  fs=paste(c('positionScoreGermline','positionScoreSomatic',
             paste('positionScoreOver',cutoff,'Somatic',sep=''),
             paste('positionScoreOver',cutoff,'Germline',sep='')),fname,sep='')

  return(fs)
  
  
}

#'Create a heatmap of tumor depth for a single-gene data frame
#'@param df is a data frame containing a single gene derived from
#'getMutationStatsForGeen
#'@param fname is name of output file
heatmapFromMutDf<-function(df=getMutationStatsForGene(gene='NF1'),fname='NF1mutations.png'){
  require(reshape2)
  require(dplyr)
  require(pheatmap)
  
  ##first start with matrix of counts
  countmat=df %>% group_by(Sample_ID,Protein_Change) %>% summarize(Count=n()) %>% acast(Sample_ID~Protein_Change,value.var='Count',fill=0)
  mut_ann=sapply(colnames(countmat),function(x) as.character(df[match(x,as.character(df$Protein_Change)),'Mutation_Type']))
  if(is.na(nrow(countmat)))
    countmat<-countmat[order(rowSums(countmat)),]
  if(is.na(ncol(countmat)))
    countmat<-countmat[,order(colSums(countmat))]
  print(dim(countmat))
    ##now separate out countmap by germline and somatic
 # par(mfrow=c(2,1))
  soms=grep('tissue',rownames(countmat))
  if(length(soms)>0){
    try(pheatmap(countmat[soms,],cluster_rows=F,cluster_cols=F,cellheight=10,cellwidth=10,annotation_col=data.frame(MutationClass=mut_ann),filename=paste('positionSomaticCount',fname,sep='')))
    try(pheatmap(countmat[-soms,],cluster_rows=F,cluster_cols=F,cellheight=10,cellwidth=10,annotation_col=data.frame(MutationClass=mut_ann),filename=paste('positionGermlineCount',fname,sep='')))
  }else{
    try(pheatmap(countmat,cluster_rows=F,cluster_cols=F,cellheight=10,cellwidth=10,annotation_col=data.frame(MutationClass=mut_ann),filename=paste('positionGermlineCount',fname,sep='')))
  
  }
  
  ##then separate out by tumor type
  types=unique(df$Mutation_Type)
  type.nums=seq(1,length(types))
  tdf=data.frame(TumorType=types,Numeric=type.nums)
  df$Numeric_Mut=tdf$Numeric[match(df$Mutation_Type,tdf$TumorType)]
  
  typemat= acast(df,Sample_ID~Protein_Change,value.var='Numeric_Mut',fill=0)
  typemat<-typemat[order(rowSums(typemat)),order(apply(typemat,2,function(x) mean(x[which(x>0)])))]
  soms=grep('tissue',rownames(typemat))
  if(length(soms)>0){
    
  pheatmap(typemat[soms,],color = c('white',rainbow(length(types))),
           cellheight=10,cellwidth=10,cluster_rows=F,cluster_cols=F,
           legend_labels=as.character(types),legend_breaks=c(type.nums)-0.5,
           filename=paste('positionTypeSomatic',fname,sep=''))
  
  pheatmap(typemat[-soms,],color = c('white',rainbow(length(types))),
           cluster_rows=F,cluster_cols=F,
           cellheight=10,cellwidth=10,
           legend_labels=as.character(types),legend_breaks=c(type.nums)-0.5,filename=paste('positionTypeGermline',fname,sep=''))
  }else{
    try(pheatmap(typemat,color = c('white',rainbow(length(types))),
             cluster_rows=F,cluster_cols=F,
             cellheight=10,cellwidth=10,
             legend_labels=as.character(types),legend_breaks=c(type.nums)-0.5,filename=paste('positionTypeGermline',fname,sep='')))
    
  }
    

}


