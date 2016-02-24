source("../../bin/clusterRNASeqData.R")


exlist=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-24/alterWgncaParams.R',wasExecuted=TRUE),
            list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/clusterRNASeqData.R',wasExecuted=TRUE))

for(signed in c(FALSE,TRUE)){
  for(sd in c(FALSE,TRUE)){
   
    
    fc.clust=clusterData(t(fc.matrix),'featureCounts')
    fc.enrich=getEnrichment(t(fc.matrix),fc.clust$tomStatic,fc.clust$TOMprefix)
    fc.eigen=evalEigenModules(t(fc.matrix),colorh1=fc.clust$tomStatic,pids=fc.pids,prefix=fc.clust$TOMprefix)
    
    
    ##then get clusters for cufflinks
    cl.clust=clusterData(t(cl.matrix),'cuffLinks')
    ##now plot eigen genes for each
    cl.enrich=getEnrichment(t(cl.matrix),cl.clust$tomStatic,cl.clust$TOMprefix)
    
    cl.eigen=evalEigenModules(t(cl.matrix),colorh1=cl.clust$tomStatic,pids=cl.pids,prefix=cl.clust$TOMprefix)
    
  }
}

#now upload feature counts files
#fcfiles=c('featureCountsWGCNAClustering.pdf','WGCNA_featureCountsTOMClusterAssignment.tsv','featureCountsTOMtop10GOTermsPermodule.csv')
#for(f in fcfiles){
#    synStore(File(f,parentId='syn5669860'), activityName='WGCNA Analysis',
#             used=c(exlist,list(list(entity='syn5051784',wasExecuted=F))))
#}




#clfiles=c('cuffLinksWGCNAClustering.pdf','WGCNA_cuffLinksTOMClusterAssignment.tsv','cuffLinksTOMtop10GOTermsPermodule.csv')
#for(f in clfiles){
#    synStore(File(f,parentId='syn5669860'),
#              used=c(exlist,list(list(entity='syn5579598',wasExecuted=F))))

#}

                                        #store on synapse
