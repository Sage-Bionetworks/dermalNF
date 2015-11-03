###
## dermalNFData.R
## This file is designed to be a basic library file so that a user can collect a specific type of dermalNF data
## Author: Sara Gosline
## Email: sara.gosline@sagebase.org
###

##SNP data files
snpfiles=synapseQuery('SELECT id,name,patientID,tissueType,tissueID FROM entity where parentId=="syn5004874"')

snpfiles<-snpfiles[grep("Final.csv",snpfiles$entity.name),]
snp.sample.names<-sapply(snpfiles$entity.name,function(x) gsub('_Final.csv','',unlist(strsplit(x,split='-'))[3]))
snp.patients<-snpfiles$entity.patientID
names(snp.patients)<-snp.sample.names

snp.tissue<-snpfiles$entity.tissueID
names(snp.tissue)<-snp.sample.names

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

    print('Now retreiving original CNV data from Dermal NF OMNI arrays...')



    sample.data<-lapply(snpfiles$entity.id,function(synid){
        print(paste("Getting sample",snpfiles$entity.name[match(synid,snpfiles$entity.id)]))
        fname=synGet(synid)
        data <- as.data.frame(fread(fname@filePath,sep=",",header=T))
        ad<-data[match(annot$Name,data$'SNP.Name'),]
        return(ad)
    })
    names(sample.data)<-snp.sample.names
    return(sample.data)
}

cnv_segmented<-function(filterSD=TRUE){

}
