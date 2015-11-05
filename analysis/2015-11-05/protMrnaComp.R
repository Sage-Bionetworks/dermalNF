#now shore up proteomics annotations, and compare with mrna


source("../../bin/dermalNFData.R")

##first assemble and store RNA-Seq data
mrna<-rna_count_files()
rna_annots<-rna_annotations()

#then collect the protin data
prot<-prot_normalized()
prot_annots<-prot_annotations()
prot_patients<-sapply(prot_annots[match(rownames(prot),prot_annots[,4]),1],gsub,'0','',fixed=T)
prot_sample<-prot_annots[match(rownames(prot),prot_annots[,4]),3]

names(prot_patients)<-rownames(prot)
names(prot_sample)<-rownames(prot)

library(pheatmap)
pheatmap(log2(prot+0.0001),annotation_row=data.frame(Patient=prot_patients,Experiment=prot_sample),cellwidth = 10,cellheight=10,file='log2ProteomicsByPatient.png')


##now reduce rna matrix to get overlap
rna.sub<-tab[rownames(tab)[which(rownames(tab)%in%colnames(prot))],]

#CHECK: why are we missing ALB and KRT1 in RNA data? They are not annotateD?  CHECK READS
zv<-which(apply(rna.sub,1,function(x) all(x==0)))
rna.sub=rna.sub[-zv,]
rna_patients<-sapply(rna_annots[match(colnames(rna.sub),rna_annots[,3]),4],function(x) gsub('0','',x))
names(rna_patients)<-colnames(rna.sub)

pheatmap(log10(rna.sub+0.001),scale='row',annotation_col=data.frame(Patient=rna_patients),cellwidth=10,cellheight=10,file='log10mRNALevelsForMeasProt.png')