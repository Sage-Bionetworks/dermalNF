##mutational data analysis


library(synapseClient)
library(data.table)
library(ggplot2)
synapseLogin()

mafdir='../../data/somatic_mafs_cleaned'
maffiles=list.files(mafdir)

toPatientId<-function(vec){
  #takes vector of tumor ids and reduces to patient
  sapply(as.character(vec),function(x) paste(unlist(strsplit(x,split='-'))[1:4],collapse='-'))
}


if(!exists('combinedMaf')){
  df1.2.maf<-gsub('.gz','',synGet('syn4924181')@filePath)
  combinedMaf<-as.data.frame(fread(df1.2.maf))
  combinedMaf$Patient<-toPatientId(combinedMaf$Tumor_Sample_Barcode)
}

if(!exists('pat.dis')){
  pat.dis<<-synTableQuery('SELECT acronym,samples from syn3281840')
  tumsByDis<<-sapply(unique(pat.dis@values$acronym),function(x) pat.dis@values$samples[which(pat.dis@values$acronym==x)])
  tumsByDis$COADREAD=c(tumsByDis$COAD,tumsByDis$READ)
}


names(maffiles)<-sapply(maffiles,function(x) unlist(strsplit(x,split='_'))[1])

getTumorIdsFromMaf<-function(maffile){
                                        #get a unique set of patients
    disname=unlist(strsplit(basename(maffile),split='_'))[1]
    tab<-as.data.frame(fread(maffile))
    pats<-toPatientId(as.character(unique(tab$Tumor_Sample_Barcode)))
    return(pats)

}

getGeneStatusByDisease<-function(disname,gene='NF1'){
    dpats<-unique(toPatientId(tumsByDis[[toupper(disname)]]))
  #  dpats<-dpats[which(names(dpats)!='PANCAN')]
    print(paste('Found',length(dpats),'samples for',disname))
    maftab<-combinedMaf[which(combinedMaf$Patient%in%dpats),]
    genetab=maftab[which(maftab$Hugo_Symbol==gene),]
    print(paste('Found',nrow(maftab),'entries for',disname,'and',nrow(genetab),'of those for',gene))
    p<-ggplot(genetab)+
      geom_bar(aes(Variant_Classification))
    return(genetab)
}

getGeneStatusFromMaf<-function(maffile,gene='NF1'){
  disname=unlist(strsplit(basename(maffile),split='_'))[1]
  tab<-as.data.frame(fread(maffile))
  genetab=tab[which(tab$gene_name==gene),]
  print(paste("Found",nrow(genetab),'entries with',gene,'mutated in',disname))
  p<-ggplot(genetab)+
  geom_bar(aes(Variant_Classification))
  return(genetab)
}

getNF1Data<-function(plot=FALSE){
    getMutDataForGene("NF1",plot)
}

getMutDataForGene<-function(geneName,plot=FALSE){
   #allnf1<-lapply(maffiles[-15],function(x) getGeneStatusFromMaf(paste(mafdir,x,sep='/'),geneName))
  allnf1<-lapply(setdiff(names(tumsByDis),'PANCAN'),function(x) getGeneStatusByDisease(x,geneName))
  names(allnf1)<-setdiff(names(tumsByDis),'PANCAN')
  dis=list()
  variant=list()
  score=list()
    tumor_id=list()
    aa_change=list()
  for(i in names(allnf1)){
    gt=allnf1[[i]]
    if(nrow(gt)==0)
      next()
    dis=c(dis,rep(i,nrow(gt)))
    variant=c(variant,gt$Variant_Classification)
    score=c(score,gt$PolyPhen)
    tumor_id=c(tumor_id,gt$Tumor_Sample_Barcode)
    aa_change=c(aa_change,gt$HGVSp)
  }

    df<-data.frame(Disease=unlist(dis),VariantClassification=unlist(variant),Score=unlist(score),Tumor=unlist(tumor_id),AAChange=unlist(aa_change))
  if(plot){
    png(paste(geneName,'mutationstatus.png',sep=''))
    p<-ggplot(df)+geom_bar(aes(VariantClassification,fill=Disease)) + ggtitle(paste(geneName,"mutations in TCGA")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
    dev.off()
  }
  return(df)
}
