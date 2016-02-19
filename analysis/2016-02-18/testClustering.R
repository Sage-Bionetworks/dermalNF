source("../../bin/clusterRNASeqData.R")

##first get clusters for feature counts
fc.clust=clusterData(t(fc.matrix),'featureCounts')


tab<-as.data.frame(fc.clust)
tab$Gene=rownames(fc.matrix)
write.table(tab,file='WGCNA_featureCountsClusterAssignment.tsv',sep='\t',row.names=F)

exlist=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-18/testClustering.R',wasExecuted=TRUE),
             list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/clusterRNASeqData.R',wasExecuted=TRUE))


fc.enrich=getEnrichment(t(fc.matrix),fc.clust$origStatic,'featureCounts')


fc.eigen=evalEigenModules(t(fc.matrix),colorh1=fc.clust$origStatic,pids=fc.pids,prefix='featureCounts')

#now upload feature counts files
fcfiles=c('featureCountsWGCNAClustering.pdf','WGCNA_featureCountsClusterAssignment.tsv','featureCountstop10GOTermsPermodule.csv')
for(f in fcfiles){
    synStore(File(f,parentId='syn5669860'), activityName='WGCNA Analysis',
             used=c(exlist,list(list(entity='syn5051784',wasExecuted=F))))
}




##then get clusters for cufflinks
cl.clust=clusterData(t(cl.matrix),'cuffLinks')

##now plot eigen genes for each
cl.enrich=getEnrichment(t(cl.matrix),cl.clust$origStatic,'cuffLinks')

cl.eigen=evalEigenModules(t(cl.matrix),colorh1=cl.clust$origStatic,pids=cl.pids,prefix='cuffLinks')


clfiles=c('cuffLinksWGCNAClustering.pdf','WGCNA_cuffLinksClusterAssignment.tsv','cuffLinkstop10GOTermsPermodule.csv')
for(f in clfiles){
    synStore(File(f,parentId='syn5669860'),
              used=c(exlist,list(list(entity='syn5579598',wasExecuted=F))))

}

                                        #store on synapse
