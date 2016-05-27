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
  all.gene.muts<<-read.table(synGet('syn6047991')@filePath,header=T,sep='\t')

#if(!exists('cancer.gene.muts')){
# cancer.genes=unique(read.table('../../data/Census_allTue Jan 19 18-58-56 2016.csv',sep=',',header=T,quote='"')$Gene.Symbol)
# cancer.gene.muts<-subset(all.gene.muts,Gene%in%cancer.genes)
#}
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


#'divideMAFfiles
#'Collects MAF files from synapse, separates them by vardict mutation status
#'@param effect is list of allowable effects.
#'@return list of list of tables
divideMAFfiles<-function(effect=c("LOW","MODERATE","MODIFIER","HIGH"),pvalthresh=0.05){

  ##goal is to filter all by effect, but only somatic/LOH by pvalue...
    allmafs<-synapseQuery("select * from entity where parentId=='syn6022474'")
    require(R.utils)
    allMuts<-lapply(allmafs$entity.id,function(x){
        res<-synGet(x)@filePath
        base=basename(res)
        f=gsub('.gz','',base)
       if(!file.exists(f)){
        file.copy(res,base)
        f=gunzip(base)[1]
        }
        tab<-fread(f,sep='\t')
#        tab<-as.data.frame(fread(input=paste('zcat < ',res@filePath)))
        etab<-subset(tab,Effect%in%effect)
        ttab<-subset(etab,`Paired-p_value`<pvalthresh)
        germ<-unique(subset(etab,Status=='Germline'))
        ss<-unique(subset(ttab,Status=='StrongSomatic'))
        ls<-unique(subset(ttab,Status=='LikelySomatic'))
        lloh<-unique(subset(ttab,Status=='LikelyLOH'))
        loh<-unique(subset(ttab,Status=='StrongLOH'))
        del<-unique(subset(ttab,Status=='Deletion'))
        print(paste("Returning tables with",nrow(germ),'germline events,',nrow(ss)+nrow(ls),'somatic events',nrow(loh)+nrow(lloh),'loh events and',nrow(del),'deletion events'))
        return(list(Germline=germ,LOH=rbind(lloh,loh),Somatic=rbind(ls,ss),Deletion=del,synId=x))
    })
    ids<-sapply(allmafs$entity.name,function(x) unlist(strsplit(x,split='.',fixed=TRUE))[1])
    names(allMuts)<-ids
    return(allMuts)
}



storeMutsForAllGenes<-function(impact=c("HIGH"),pvalthresh=0.05){
    mafs=divideMAFfiles(impact,pvalthresh=pvalthresh)
    synids<-lapply(mafs,function(x) x$synId)
    all.mafs<-do.call('rbind',lapply(mafs,function(x) rbind(x$Germline,x$Somatic,x$LOH,x$Deletion)))
  fname=paste('varDictMutations_AllGenesAllSamples_pval',pvalthresh,'_impact',paste(impact,collapse='_'),'.tsv',sep='')
  write.table(all.mafs,file=fname, sep='\t',row.names=F,quote=F,col.names=T)
  #gfname<-paste(fname,'.gz')
                                        #gzip(fname,gfname)
      ul<-lapply(as.character(synids),function(x) list(entity=x))
#  ul<-lapply(synids,function(x) list(entity=x))
  synStore(File(fname,parentId='syn5605256'),used=ul,executed='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/WGSData_VarDict.R')

}


#'getMutationStatsForGene obtains all mutations for a particular gene of interest across all patients
#'@param gene is the gene symbol in question
#'@param impact is a list of which mutations to include, defaults to all ('HIGH','MODERATE' and 'LOW')
#'@param doPlot: if set to true, will plot some basic statistics about where and when this mutation occurs
getMutationStatsForGene<-function(expr.gene.muts,gene='NF1',doPlot=FALSE,effect=c("LOW","MODERATE","HIGH"),prefix='p05'){
  
  sdf<-subset(expr.gene.muts,Gene==gene)
  # if()
  sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  sdf$Patient<-sapply(sdf$Patient,function(x) paste(x,'(n=',length(unique(sdf$Sample[which(sdf$Patient==x)])),')'))
  transcripts<-unique(sdf$Transcript)
  print(paste("Found",length(transcripts),'unique transcripts for',gene))
  sdf$Effect=factor(as.character(sdf$Effect),levels=c("LOW","MODERATE",'HIGH'))
  sdf<-subset(sdf,Effect%in%effect)
  if(nrow(sdf)==0)
    return(sdf)
  
  if(doPlot){
    for(t1 in transcripts){
      tdf<-subset(sdf,Transcript==t1)
      p<-ggplot(tdf)+geom_jitter(aes(y=Patient,x=Codon_Change,col=Functional_Class,shape=Status,size=Effect))
      #   p  +facet_grid(.~Status)
      p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(paste('Gene',gene,'Transcript',t1))#+facet_grid(.~Status)
      ggsave(p,file=paste('gene',gene,'transcript',t1,'mutationsByType',prefix,'.png',sep=''))
    }
  }
  return(sdf) 
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
