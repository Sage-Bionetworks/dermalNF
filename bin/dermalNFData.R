###
## dermalNFData.R
## This file is designed to be a basic library file so that a user can collect a specific type of dermalNF data
## Author: Sara Gosline
## Email: sara.gosline@sagebase.org
###
library(synapseClient)
synapseLogin()
library(data.table)
require(parallel)

#################
# CNV
#################

cnv_annotations<-function(){
    snpfiles=synapseQuery('SELECT id,name,patientID,tissueType,tissueID,alternateTumorID FROM entity where parentId=="syn5004874"')

    names(snpfiles)<-c('tissueType','patientId','alternateTumorId','File','tissueId','synapseId')
    return(snpfiles)
}

cnv.dat<-cnv_annotations()
cnv.dat<-cnv.dat[which(!is.na(cnv.dat$patientId)),]
patients<-cnv.dat$patientId
names(patients)<-sapply(cnv.dat$File,function(x) gsub('3096-PBK-','X',gsub('_Final.csv','',x)))
patients<-sapply(patients,function(x) gsub("CT0*","",x))
names(patients)<-cnv.dat$synapseId
clnames<-paste(patients,cnv.dat$tissueId)
names(clnames)<-names(patients)

tissueType=cnv.dat$tissueType
names(tissueType)<-names(patients)


#SNP annotation file
snp_annotation_file<-function(){
  ##need to downlod and read in large annotation file as well
  print("Retrieving OMNI Array SNP annotation data from Synapse...")
  anndata<-synGet('syn5297573')
  return(anndata@filePath)
}

snp_annotation_data<-function(){
  fp=snp_annotation_file()
    annot <- as.data.frame(fread(fp,sep=",",header=T))
    return(annot)
}

cnv_unprocessed_files<-function(){

  snpfiles=synapseQuery('SELECT id,name,patientID,tissueType,tissueID FROM entity where parentId=="syn5004874"')
  snpfiles<-snpfiles[grep("Final.csv",snpfiles$entity.name),]
  snp.sample.names<-sapply(snpfiles$entity.name,function(x) gsub('_Final.csv','',unlist(strsplit(x,split='-'))[3]))
  snp.patients<-snpfiles$entity.patientID
  names(snp.patients)<-snp.sample.names

  snp.tissue<-snpfiles$entity.tissueID
  names(snp.tissue)<-snp.sample.names

    sample.data<-lapply(snpfiles$entity.id,function(synid){
    print(paste("Getting sample",snpfiles$entity.name[match(synid,snpfiles$entity.id)]))
    fname=synGet(synid)
    return(fname@filePath)
  })
    names(sample.data)<-snpfiles$entity.id
    return(sample.data)
}

#this function gets the original files from the OMNI arrays
cnv_unprocessed<-function(annot=NA){
    if(is.na(annot))
        annot=snp_annotation_data()

    ##SNP data files
    snpfiles=synapseQuery('SELECT id,name,patientID,tissueType,tissueID FROM entity where parentId=="syn5004874"')

    snpfiles<-snpfiles[grep("Final.csv",snpfiles$entity.name),]
    snp.sample.names<-sapply(snpfiles$entity.name,function(x) gsub('_Final.csv','',unlist(strsplit(x,split='-'))[3]))
    snp.patients<-snpfiles$entity.patientID
    names(snp.patients)<-snp.sample.names

    snp.tissue<-snpfiles$entity.tissueID
    names(snp.tissue)<-snp.sample.names


    print('Now retreiving original CNV data from Dermal NF OMNI arrays...')
    #here get the sample data from snp files
    sample.data<-lapply(snpfiles$entity.id,function(synid){
        print(paste("Getting sample",snpfiles$entity.name[match(synid,snpfiles$entity.id)]))
        fname=synGet(synid)
        data <- as.data.frame(fread(fname@filePath,sep=",",header=T))
        ad<-data[match(annot$Name,data$'SNP.Name'),]
        return(ad)
    })
    names(sample.data)<-snpfiles$entity.id
    return(sample.data)
}

#this gets the CNV segment files
cnv_segmented<-function(filterSD=TRUE){
    if(filterSD)
        si='syn5049753'
    else
        si='syn5049755'
    fn<-synGet(si)
    tab<-read.table(fn@filePath,header=T)
    return(tab)

}

cnv_segmented_by_gene<-function(){
    si='syn5462050'
    fn<-synGet(si)
    tab<-read.table(fn@filePath,header=T)
    return(tab)
}

#################
#PROTEOMICS
#################
protein_annotations<-function(){
    annots<-synapseQuery("select name,ID,dataType,tissueID,tissueType,patientID,sampleID from entity where parentId=='syn4984949'")
    annots<-annots[-grep('EMPTY',annots$entity.name),]

    colnames(annots)<-c('tissueType','dataType','sampleId','patientId','fileName','tissueId','synapseId')

    return(annots)
}

#this merely calculates the ratio for each file
get.protein.from.file<-function(sn,top_only=FALSE){
    sd<-synGet(unlist(sn))
    tab<-read.table(sd@filePath,sep='\t',header=T,as.is=T)
    nums<-tab[,6]
    denoms<-tab[,7]
    ratios<-nums/denoms

  groups<-tab[,1]

  u.groups<-unique(groups)
  groups.ids<-sapply(u.groups,function(x) return(paste(unique(tab[which(tab[,1]==x),4]),collapse=';')))
  names(groups.ids)<-u.groups

  print(paste("Found",length(u.groups),'unique protein groups'))

  u.tops<-sapply(u.groups,function(x) intersect(which(tab[,1]==x),which(tab[,3]=='TOP PROTEIN')))
  #now filter for top
  top.ratios=ratios[u.tops]
  top.nums=nums[u.tops]
  top.conts=denoms[u.tops]

  return(list(Ratios=top.ratios,Raw=top.nums,Control=top.conts,Prot.ids=tab[u.tops,4],Origin=tab[,8]))

}


prot_normalized<-function(store=FALSE,all.expr=TRUE){
  #store indicates we should calculate the values and uplod to synapse, otherwise we can just download pre-computed
  #all.expr means select only those proteins that non-zero in at least one sample

  allfiles= synapseQuery('SELECT name,ID,patientID,tissueID FROM entity WHERE parentId=="syn4984949"')



  if(store){
    res<-sapply(allfiles$entity.id,function(x) get.protein.from.file(x,TRUE))
    names(res)<-allfiles$entity.id
    #first col  lect all proteins annotated in any file
    all.prots<-NULL
    for(i in 1:ncol(res))
        all.prots<-union(all.prots,res[['Prot.ids',i]])
    #filter for those that are expressed across all samples
 #   expr.prots<-res[['Prot.ids',1]]
#    for(i in 2:ncol(res))
#        expr.prots<-intersect(expr.prots,res[['Prot.ids',i]])

    prot.ids<-unique(unlist(sapply(all.prots,function(x) unlist(strsplit(x,split=';')))))

                                        #now create biomart mapping
    require(biomaRt)
    ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host='www.ensembl.org')
    filters = listFilters(ensembl)
    attributes = listAttributes(ensembl)

    epep="ensembl_peptide_id"
    egene='hgnc_symbol'
    gene.mapping<-getBM(attributes=c(epep,egene),filters=c(epep),values=as.list(prot.ids),mart=ensembl)

    allsamps<-colnames(res)
    expr.ratio.mat<-sapply(all.prots,function(x){
       # pvec<-NULL
                                        # samps<-NULL

        pvec<-sapply(allsamps,function(i){
            rv<-grep(x,res[['Prot.ids',i]])
            if(length(rv)==0)
                return(0)
            else
                return(res[['Ratios',i]][rv])
        })
        names(pvec)<-allsamps
        unlist(pvec)
    })

    gn<-gene.mapping[match(colnames(expr.ratio.mat),gene.mapping[,1]),2]
    expr.ratio.mat[which(is.na(expr.ratio.mat),arr.ind=T)]<-0.0
    #expr.ratio.mat<-expr.ratio.mat[-grep('EMPTY',rownames(expr.ratio.mat)),]
    gn<-gene.mapping[match(colnames(expr.ratio.mat),gene.mapping[,1]),2]
    gn[which(is.na(gn))]<-colnames(expr.ratio.mat)[which(is.na(gn))]
    colnames(expr.ratio.mat)<-gn
    expr.ratio.mat<-t(expr.ratio.mat)
    df=data.frame(Protein=rownames(expr.ratio.mat),expr.ratio.mat)
    write.table(df,file='proteinFoldChangeOverControl.txt',sep='\t',row.names=F)
    print('Storing file on Synapse...')
    synStore(File('proteinFoldChangeOverControl.txt',parentId='syn4984703'),
             used=list(c(sapply(allfiles$entity.id,function(x) list(entity=x)),list(name='dermalNFData.R',url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/dermalNFData.R'))),
             activityName='Computed ratios between protein and control',
             activityDescription='called prot_normalized with store=TRUE')
    expr.ratio.mat<-df
  }else{
    matfile=synGet('syn5305003')
    expr.ratio.mat<-as.data.frame(fread(matfile@filePath,sep='\t',header=T))
  }

  if(all.expr){
    zo=which(apply(expr.ratio.mat[,-1],1,function(x) all(x==0)))
    if(length(zo)>0){
      print(paste('Removing',length(zo),'proteins from matrix because they have only 0-values'))
      expr.ratio.mat<-expr.ratio.mat[-zo,]
    }
    return(expr.ratio.mat)
  }

  return(expr.ratio.mat)
}



#################
#RNA
#################
rna_annotations<-function(){
    synq=synapseQuery("select name,id,patientID,tissueID,alternateTumorID from entity where parentId=='syn5493036'")
    colnames(synq)<-c('patientId','alternateTumorId','fileName','tissueId','synapseId')
    synq=synq[grep('_featureCounts.txt',synq$fileName),]
    return(synq)
}

rna_bam_annotations<-function(){
    synq=synapseQuery("select name,id,patientID,tissueID,alternateTumorID from entity where parentId=='syn4984620'")
    colnames(synq)<-c('patientID','alternateTumorId','fileName','tissueId','synapseId')
    synq=synq[grep('.bam$',synq$fileName),]
    return(synq)
}

##here are the count files analyzed by featureCounts
rna_count_matrix<-function(stored=TRUE,doNorm=FALSE,minCount=0,doLogNorm=FALSE){

    if(!stored){
        synq=synapseQuery("select name,id,patientID,tissueID from entity where parentId=='syn5493036'")
        synq<-synq[grep("accepted_hits",synq$entity.name),]
        synfiles<-sapply(synq$entity.id,synGet)
                                        #now read in alfilel values

        allfs<-lapply(synfiles,function(x) read.table(x@filePath,header=T,as.is=T))
        names(allfs)<-synq$entity.id

                                        #now get individual genes to create data matrix
        hugo.genes<-unique(allfs[[1]][,2])


                                        #now let's get individual counts across patient samples
        gene.pat.mat<-sapply(hugo.genes,function(x,allfs){
            res<-sapply(names(allfs),function(y){
                mat<-allfs[[y]]
                sum(mat[which(mat[,2]==x),1])})
            names(res)<-names(allfs)
            res
        },allfs)

        colnames(gene.pat.mat)<-hugo.genes

        write.table(gene.pat.mat,file='featureCountsByGeneBySample.txt',row.names=T,col.names=T)
        sf=File('featureCountsByGeneBySample.txt',parentId='syn4984701')
        synStore(sf,used=list(list(name='dermalNFData.R',
                        url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/dermalNFData.R')),
                 activityName='Create matrix of all counts across samples')
    }else{
        gene.pat.mat<-read.table(synGet('syn5051784')@filePath)
    }

    gene.pat.mat<-t(gene.pat.mat)

    if(doNorm){
      print('Performing size factor adjustment to samples')
      require(DESeq2)
      samp=data.frame(SampleID=colnames(gene.pat.mat))
      cds<- DESeqDataSetFromMatrix(gene.pat.mat,colData=samp,~SampleID)#now collect proteomics data

      sizeFac<-estimateSizeFactors(cds)

      normCounts<-assay(cds)/sizeFac@colData$sizeFactor
      colnames(normCounts)<-colnames(gene.pat.mat)
      gene.pat.mat<-normCounts

    }else if(doLogNorm){
      print("Performing variance stabilizing log2 normalization")
      require(DESeq2)
      samp=data.frame(SampleID=colnames(gene.pat.mat))
      cds<- DESeqDataSetFromMatrix(gene.pat.mat,colData=samp,~SampleID)#now collect proteomics data
      vstab=rlog(cds)

      varmat<-assay(vstab)
      colnames(varmat)<-colnames(gene.pat.mat)
      gene.pat.mat<-varmat
      minCount=log2(minCount)

    }

    sel.vals=which(apply(gene.pat.mat,1,function(x) all(x>=minCount)))

    return(gene.pat.mat[sel.vals,])

}


#we can also get the FPKM
rna_fpkm_matrix<-function(byIsoform=FALSE){
  ##DOES NOT WORK YET....
    if(byIsoform){
        gene.pat.mat<-read.table(synGet('syn5579597')@filePath)
    }
    else{
        gene.pat.mat<-read.table(synGet('syn5579598')@filePath,row.names=NULL)
        dupes<-unique(gene.pat.mat[which(duplicated(gene.pat.mat[,1])),1])
        dupe.vals<-t(sapply(dupes,function(x)
            colSums(gene.pat.mat[which(gene.pat.mat[,1]==x),2:ncol(gene.pat.mat)])))
        sing.vals<-gene.pat.mat[which(!gene.pat.mat[,1]%in%dupes),2:ncol(gene.pat.mat)]
        rownames(dupe.vals)<-dupes
        rownames(sing.vals)<-gene.pat.mat[which(!gene.pat.mat[,1]%in%dupes),1]
        newdf<-rbind(dupe.vals,sing.vals)
        gene.pat.mat<-newdf
    }
    #gene.pat.mat<-t(gene.pat.mat)
    return(gene.pat.mat)
}

#################
#WGS
#################
wgs_annotations<-function(){
    synq=synapseQuery("select name,id,patientID,tissueID,alternateTumorID from entity where parentId=='syn5522788'")
    colnames(synq)<-c('patientId','alternateTumorId','fileName','tissueId','synapseId')
   # synq=synq[-grep('hard-filtered',synq$fileName),]
    return(synq)
}
