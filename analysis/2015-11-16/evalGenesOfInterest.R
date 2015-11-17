##here we collect a list of genes and evaluate their expression in the proteomics and RNA samples

genelist<-c('RET','GFRA1','GDNF','NRTN','ARTN','PSPN','NGFR',
            'GAP43','NTRK1','NTRK2','NTRK3','NGF','BDNF','NTF3','NTF4',
            'UCHL1','S100B','ENO2')

source("../../bin/dermalNFData.R")

#get rna data
rna.annot<-rna_annotations()
rna<-rna_count_matrix(doNorm=TRUE,minCount=5)

allvars=apply(rna,1,function(x) var(x,na.rm=T)^2/mean(x,na.rm=T))
pats<-rna.annot$patientId
names(pats)<-rna.annot$synapseId

library(pheatmap)
pheatmap(log10(0.00001+rna[order(allvars,decreasing=T)[1:100],]),annotation_col=data.frame(Patient=pats),
         cellwidth=10,cellheight=10,
          clustering_distance_cols = 'correlation',clustering_distance_rows = 'correlation',
         file='RNA_100_mostVariable_min5counts.png')

colnames(rna)<-paste(pats,rna.annot$tissueID)
plot(hclust(dist(t(rna)),method='ward.D2'))


#now get all mrnAs, even the poorly expressed
rna<-rna_count_matrix(doNorm=TRUE,minCount=0)

red.set<-rna[match(genelist,rownames(rna)),]
pheatmap(log10(0.00001+red.set),annotation_col=data.frame(Patient=pats),
         clustering_distance_cols = 'correlation',file='signalingProteinGeneExpress.png')

#now get the protein data
prot<-prot_normalized()
prot.annot=protein_annotations()
pats=gsub('CT0*','',prot.annot$patientId)
names(pats)<-prot.annot$synapseId
exper=gsub('dermalNF_proteomics_sample_X','run_',gsub('_Normalized_P_SC__NP_SC_.txt','',prot.annot$fileName))
names(exper)<-prot.annot$synapseId
pheatmap(log10(prot[,-1]+0.0001),annotation_col=data.frame(Patients=pats),
         cellheight = 10,cellwidth = 10,labels_row=prot[,1],file='log10_nonzero_Proteins.png',
         clustering_distance_rows='correlation',clustering_distance_cols='correlation',
         clustering_method='ward.D2')

synStore(File('log10_nonzero_Proteins.png',parentId='syn4984701'),
         activityName='Log10 heatmap of proteins that are expressed in at least one experiment',
         used=list(list(name='evalGenesOfInterest.R',url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2015-11-16/evalGenesOfInterest.R')))
      
 # pheatmap(log10(prot[,-1]+0.0001),annotation_col=data.frame(Patients=pats),cellheight = 10,cellwidth = 10,labels_row=prot[,1],file='allProteins.png')
##