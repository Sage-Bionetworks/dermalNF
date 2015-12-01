###compare RNAseq data to CNV/proteomics data
##Goal of this is to create universal mapping across all the different types of data.
##TODO:  include mapping for analysis files as well as data files

source('../../bin/dermalNFData.R')







getPatientMapping<-function(patientId=NA){
    #here we can get the patient mappings from the synapse talbe
    schemaId='syn5478651'
    if(is.na(patientId))
        queryResult<-synTableQuery(sprintf("select * from %s", schemaId))
    else
        queryResult<-synTableQuery(sprintf("select * from %s where patient='%s'", schemaId,patientId))
    return(queryResult@values)

}

storePatientMappingDataInTable<-function(){
                                        #this function is the primary one that performs the patient sample mapping to add into
                                        #a synapse table

                                        #first get all the annotations
    rna.map<-rna_bam_annotations()
    cnv.map<-cnv_annotations()
    wgs.map<-wgs_annotations()
    prot.map<-protein_annotations()

    ##next we need to get the mappings from one sample to another.
    ##this will be in the annotations...
    ##for each patient, for each sample, identify cnv, proteomics, RNA and WGS samples
    ##TODO: map the analysis files as well.

    all.patients<-as.character(seq(1,13))
    sample.list<-NULL
                                        #patient.sample.matches<-sapply(all.patients,function(x){\
    for(x in all.patients){
        #TODO: harmonize identifiers across all...
        rna.samps<-unique(rna.map[grep(paste("CT0+",x,'$',sep=''),rna.map$patientID),])
        dna.samps<-unique(wgs.map[grep(paste("CT0+",x,'$',sep=''),wgs.map$patientID),])
        prot.samps<-unique(prot.map[grep(paste("CT0+",x,'$',sep=''),prot.map$patientId),])
        cnv.samps<-unique(cnv.map[grep(paste("CT0+",x,'$',sep=''),cnv.map$patientID),])

        ##add in all samples that match
        matched.rna<-intersect(rna.samps$tissueID,dna.samps$alternateTumorID)
        print(paste('Found',length(matched.rna),'RNA samples that match dna samples for patient',x))
        pcount=1
        for(mr in matched.rna){
            dna.tiss=dna.samps[match(mr,dna.samps$alternateTumorID),'tissueID']
            wgs.si=dna.samps[match(mr,dna.samps$alternateTumorID),'synapseId']
            prot.si<-prot.samps[match(dna.tiss,prot.samps$tissueId),'synapseId']
            cnv.si<-cnv.samps[match(dna.tiss,cnv.samps$tissueID),'synapseID']
            rna.si<-rna.samps[match(mr,rna.samps$tissueID),'synapseID']
                                        # sample.list[[paste(x,pcount,sep='_')]]<-
            sample.list<-rbind(sample.list,list(patient=x,DnaID=dna.tiss,RnaID=mr,RNASeq=rna.si,WGS=wgs.si,Prot=prot.si,CNV=cnv.si))
            pcount<-pcount+1
        }

        #next add in unmatched dna
        unmatched.dna<-setdiff(dna.samps$tissueID,rna.samps$alternateTumorID)
        for(mr in unmatched.dna){
            dna.tiss=mr
            wgs.si=dna.samps[match(dna.tiss,dna.samps$tissueID),'synapseId']
            prot.si<-prot.samps[match(dna.tiss,prot.samps$tissueId),'synapseId']
            cnv.si<-cnv.samps[match(dna.tiss,cnv.samps$tissueID),'synapseID']

            rna.si<-NA
                                        #      sample.list[[paste(x,pcount,sep='_')]]<-
            sample.list<-rbind(sample.list,list(patient=x,DnaID=dna.tiss,RnaID=NA,RNASeq=rna.si,WGS=wgs.si,Prot=prot.si,CNV=cnv.si))
            pcount<-pcount+1
        }

        ##then add in unamtched rna
        unmatched.rna<-setdiff(rna.samps$tissueID,dna.samps$alternateTumorID)
        for(mr in unmatched.rna){
            wgs.si=NA
            prot.si<-NA
            cnv.si<-NA
            rna.si<-rna.samps[which(rna.samps$tissueID==mr),'synapseID']

                                        #      sample.list[[paste(x,pcount,sep='_')]]<-
            sample.list<-rbind(sample.list,list(patient=x,DnaID=NA,RnaID=mr,RNASeq=rna.si,WGS=wgs.si,Prot=prot.si,CNV=cnv.si))
            pcount<-pcount+1
        }

    }
    #create data frame and unlist
    samp.df<-as.data.frame(sample.list)
    for(i in 1:ncol(samp.df))
        samp.df[,i]<-unlist(samp.df[,i])
                                        #now try to store as synapse table

    tres<-as.tableColumns(samp.df)
    cols<-tres$tableColumns
    fileId=tres$fileHandleId
    projectId='syn4984604'
    schema<-TableSchema(name='Dermal NF Data Summary',parent=projectId,columns=cols)
    table<-Table(schema,fileId)
    table<-synStore(table,retrieveData=TRUE)
    return(table)
}
