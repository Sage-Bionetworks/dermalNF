##use the featureCounts package to re-process the BAMs
##need to uplod this script and add in the annotations more effectively
library(Rsubread)
library(synapseClient)

synapseLogin()

synfiles=synapseQuery("select id,name,Patient_ID,Tissue_ID from entity where parentId=='syn4984620'")
#get all bam files

#removed the unmapped
bamfiles<-grep("accepted_hits_SL[0-9]*.bam$",synfiles[,2])

                                        #gene.mapping file
xref=synGet('syn5014885',downloadFile=T)
gtf=synGet('syn5014886',downloadFile=T)
xref<-read.table(xref@filePath,header=F,as.is=T,sep='\t',quote='"')
gene.mapping<-xref[,5]
names(gene.mapping)<-xref[,1]

library(parallel)

#re-quantify each file - this can take a while!
allfs<-mclapply(bamfiles,function(x,gene.mapping){
    res<-featureCounts(paste(bamfiledir,x,sep='/'),
                       annot.ext=gtf@filePath,isGTFAnnotationFile=TRUE,isPairedEnd=TRUE)
    new.tab<-data.frame(Counts=res$counts,Symbol=sapply(rownames(res$counts),function(y) gene.mapping[[y]]))
    new.tab
},gene.mapping,mc.cores=2)

names(allfs)<-sapply(bamfiles,function(x) unlist(strsplit(x,split='.',fixed=T))[1])

                                        #write to table
ucsc.genes<-unique(rownames(allfs[[1]]))
hugo.genes<-unique(allfs[[1]][,2])

sapply(names(allfs),function(x){
    write.table(allfs[[x]],file=paste(x,'_featureCounts.txt',sep=''),row.names=T,col.names=T)
})
