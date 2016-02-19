source("../../bin/clusterRNASeqData.R")

##first get clusters for feature counts
fc.clust=clusterData(t(fc.matrix),'featureCounts')


fc.enrich=getEnrichment(t(fc.matrix),fc.clust$origStatic,'featureCounts')

fc.eigen=evalEigenModules(t(fc.matrix),colorh1=fc.clust$origStatic,pids=fc.pids,prefix='featureCounts')


##then get clusters for cufflinks
cl.clust=clusterData(t(cl.matrix),'cuffLinks')

##now plot eigen genes for each
cl.enrich=getEnrichment(t(cl.matrix),cl.clust$origStatic,'featureCounts')

cl.eigen=evalEigenModules(t(cl.matrix),colorh1=cl.clust$origStatic,pids=cl.pids,prefix='cuffLinks')


##now we should do functional enrichment


