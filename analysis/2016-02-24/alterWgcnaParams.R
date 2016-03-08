source("../../bin/clusterRNASeqData.R")


exlist=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-02-24/alterWgncaParams.R',wasExecuted=TRUE),
            list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/clusterRNASeqData.R',wasExecuted=TRUE))

for(signed in c(FALSE,TRUE)){
  for(sd in c(FALSE,TRUE)){

   fc.clust=clusterData(t(fc.matrix),'featureCounts',signed,sd,topGenes=3000)
  #  fc.enrich=getEnrichment(fc.clust$expr,fc.clust$tomStatic,fc.clust$TOMprefix)
    fc.eigen=evalEigenModules(fc.clust$expr,colorh1=fc.clust$tomStatic,pids=fc.pids,prefix=fc.clust$TOMprefix)
    write.table(fc.eigen,paste('featureCounts',ifelse(sd,'sd','conn'),'filtered',ifelse(signed,'signed','unsigned'),'clusterEigenGenes.tab',sep='_'))

    ##then get clusters for cufflinks
    cl.clust=clusterData(t(cl.matrix),'cuffLinks',signed,sd,topGenes=3000)
    ##now plot eigen genes for each
    #cl.enrich=getEnrichment(cl.clust$expr,cl.clust$tomStatic,cl.clust$TOMprefix)
    cl.eigen=evalEigenModules(cl.clust$expr,colorh1=cl.clust$tomStatic,pids=cl.pids[rownames(cl.clust$expr)],prefix=cl.clust$TOMprefix)
    write.table(cl.eigen,paste('cuffLinks',ifelse(sd,'sd','conn'),'filtered',ifelse(signed,'signed','unsigned'),'clusterEigenGenes.tab',sep='_'))  }
}


all.files=list.files('.')
go.files=all.files[grep('GO',all.files)]

all.files<- all.files[grep('Cluster',all.files,ignore.case=T)]

fcfiles=c(all.files[grep('*featureCounts*',all.files)],go.files[grep('*featureCounts*',go.files)])
clfiles=c(all.files[grep('*cuffLinks*',all.files)],go.files[grep('cuffLinks',go.files)])

#now upload feature counts files
#fcfiles=c('featureCountsWGCNAClustering.pdf','WGCNA_featureCountsTOMClusterAssignment.tsv','featureCountsTOMtop10GOTermsPermodule.csv')
for(f in fcfiles){
    synStore(File(f,parentId='syn5669860'), activityName='WGCNA Analysis',
             used=c(exlist,list(list(entity='syn5051784',wasExecuted=F))))
}




#clfiles=c('cuffLinksWGCNAClustering.pdf','WGCNA_cuffLinksTOMClusterAssignment.tsv','cuffLinksTOMtop10GOTermsPermodule.csv')
for(f in clfiles){
    synStore(File(f,parentId='syn5669860'),
              used=c(exlist,list(list(entity='syn5579598',wasExecuted=F))))

}

                                        #store on synapse
