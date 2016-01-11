###compare RNAseq data to CNV/proteomics data
##this is primarily a test file for the data mapping file...
source('../../bin/dermalNFData.R')

#mapped to gene
rna_counts=rna_count_matrix()

#mapped to genes
proteomics=prot_normalized()

#so this data isn't mapped to a gene yet
#cnv=cnv_segmented()

cnv<-cnv_segmented_by_gene()

#now get all the annotations
rna.map<-rna_bam_annotations()
cnv.map<-cnv_annotations()
wgs.map<-wgs_annotations()
prot.map<-protein_annotations()

##last we need to get the mappings from one sample to another. this will be in the annotations...
##for each patient, for each sample, identify cnv, proteomics, RNA and WGS samples

all.patients<-as.character(seq(1,13))
sample.list<-NULL
#patient.sample.matches<-sapply(all.patients,function(x){\
for(x in all.patients){
    rna.samps<-unique(rna.map[grep(paste("CT0+",x,'$',sep=''),rna.map$patientID),])
    dna.samps<-unique(wgs.map[grep(paste("CT0+",x,'$',sep=''),wgs.map$patientID),])
    prot.samps<-unique(prot.map[grep(paste("CT0+",x,'$',sep=''),prot.map$patientId),])
    cnv.samps<-unique(cnv.map[grep(paste("CT0+",x,'$',sep=''),cnv.map$patientID),])
    
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

##now select patient samples