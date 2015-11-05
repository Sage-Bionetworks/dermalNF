###
## dermalNFData.R
## This file is designed to be a basic library file so that a user can collect a specific type of dermalNF data
## Author: Sara Gosline
## Email: sara.gosline@sagebase.org
###
library(synapseClient)
synapseLogin()




#################
# CNV
#################

cnv_annotation<-function(){
    snpfiles=synapseQuery('SELECT id,name,patientID,tissueType,tissueID FROM entity where parentId=="syn5004874"')

    names(snpfiles)<-c('tissueType','patientID','File','tissueID','synapseID')
    return(snpfiles)
}

#SNP annotation file
snp_annotation_data<-function(){
    ##need to downlod and read in large annotation file as well
    print("Retrieving OMNI Array SNP annotation data from Synapse...")
    anndata<-synGet('syn5005069')
    annot <- as.data.frame(fread(anndata@filePath,sep=",",header=T))
    return(annot)
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

#################
#PROTEOMICS
#################
protein_annotation<-function(){
    annots<-synapseQuery("select name,ID,dataType,tissueID,tissueType,patientID,sampleID from entity where parentID=='syn4984949'")
    colnames(annots)<-c()

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

prot_annotations<-function(){
    allfiles= synapseQuery('SELECT name,ID,patientID,tissueID FROM entity WHERE parentId=="syn4984949"')
    colnames(allfiles)<-c('')
    return(allfiles)
}

prot_normalized<-function(all.expr=TRUE){
    allfiles= synapseQuery('SELECT name,ID,patientID,tissueID FROM entity WHERE parentId=="syn4984949"')

    res<-sapply(allfiles$entity.id,function(x) get.protein.from.file(x,TRUE))
    names(res)<-allfiles$entity.id

    #first collect all proteins annotated in any file
    all.prots<-NULL
    for(i in 1:ncol(res))
        all.prots<-union(all.prots,res[['Prot.ids',i]])
    #filter for those that are expressed across all samples
    expr.prots<-res[['Prot.ids',1]]
    for(i in 2:ncol(res))
        expr.prots<-intersect(expr.prots,res[['Prot.ids',i]])


    prot.ids<-unique(unlist(sapply(all.prots,function(x) unlist(strsplit(x,split=';')))))

    if(all.expr)
        all.prots<-expr.prots
                                        #now create biomart mapping
    require(biomaRt)
    ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
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
    return(expr.ratio.mat)



}



#################
#RNA
#################
rna_annotations<-function(){
    synq=synapseQuery("select name,id,Patient_ID,Tissue_ID from entity where parentId=='syn4984701'")
    colnames(synq)<-c('tissueID','fileName','synapseId','patientId')
    return(synq)
}

rna_bam_files<-function(){

}


##here are the count files analyzed by featureCounts
rna_count_matrix<-function(stored=TRUE){

    if(!stored){
        synq=synapseQuery("select name,id,Patient_ID,Tissue_ID from entity where parentId=='syn4984701'")
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
    return(gene.pat.mat)

}


#################
#WGS
#################
