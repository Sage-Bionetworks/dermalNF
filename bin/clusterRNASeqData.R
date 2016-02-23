#'
#' Look for patterns of interest across patients from the RNA
#' Goal is to find networks of 'association' that contain tightly clustered genes and/or pathways
#'
#'
#'

source("../../bin/dermalNFData.R")

library(WGCNA)
options(stringsAsFactors = FALSE)
library(cluster)

##collect two matrices, one with feature counts and one from cufflinks
fc.matrix = rna_count_matrix(stored=TRUE,doNorm=TRUE,minCount=2,doLogNorm=FALSE)
patients=rna_annotations()[,c('patientId','synapseId')]

fc.pids=sapply(patients$patientId,function(x) gsub('CT0+','',x))
names(fc.pids)<-patients$synapseId

cl.matrix = rna_fpkm_matrix(byIsoform=FALSE)
vars=apply(cl.matrix,1,var,na.rm=T)
cl.matrix=cl.matrix[-which(vars==0),]
cl.patients=fpkm_annotations()[,c('patient','sample')]
cl.pids=cl.patients$patient
names(cl.pids)<-cl.patients$sample

#' do WGCNA analysis to get clusters of gene modules
#'
clusterData <- function(datExpr,prefix=''){

    ADJ1=abs(cor(datExpr,use="p"))^6
    k=softConnectivity(datE=datExpr,power=6)
    pdf(paste(prefix,'WGCNAClustering.pdf',sep=''),10,5)

    #sizeGrWindow(10,5)
    par(mfrow=c(1,2))
    hist(k)
    scaleFreePlot(k, main="Check scale free topology\n")

    fullDatExpr=datExpr
    datExpr=datExpr[, rank(-k,ties.method="first" )<=3600]
    dissADJ=1-ADJ1


    hierADJ=hclust(as.dist(dissADJ), method="average" )

    ###CUTTING TREEE
    ##first experiment with various ways of labeling branches
    branch.number=cutreeDynamic(hierADJ,method="tree")
    # This function transforms the branch numbers into colors
    colorDynamicADJ=labels2colors(branch.number )

   #sizeGrWindow(10,5);
  #  plotDendroAndColors(hierADJ, colors=colorDynamicADJ, dendroLabels = FALSE, hang = 0.03,
   #                     main = "Gene hierarchical clustering dendrogram and assigned colors" )

    ##now restrict by module size
    colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
    # Plot the dendrogram with module colors
   # sizeGrWindow(10,5);
  #  plotDendroAndColors(hierADJ, colors = data.frame(colorStaticADJ),
  #                      dendroLabels = FALSE, abHeight = 0.99,
  #                      main = "Gene dendrogram and module colors")


    ##now do the hybrid approach
    colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ,
                                                      cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))
   #  sizeGrWindow(10,5)
    plotDendroAndColors(dendro = hierADJ,
                        colors=data.frame(colorStaticADJ,
                                          colorDynamicADJ, colorDynamicHybridADJ),
                        dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                        main = "Gene dendrogram and module colors")

    ###TOM analysis
    #topology based distance
    dissTOM=TOMdist(ADJ1)
    hierTOM = hclust(as.dist(dissTOM),method="average")

    colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
    colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
    colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                        deepSplit=2, pamRespectsDendro = FALSE))
    # Now we plot the results
   # sizeGrWindow(10,5)
    plotDendroAndColors(hierTOM,
                        colors=data.frame(colorStaticTOM,
                                          colorDynamicTOM, colorDynamicHybridTOM),
                        dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                        main = "Gene dendrogram and module colors, TOM dissimilarity")

    dev.off()
    ##which clustering should i return?
    ret=list(origStatic=colorStaticADJ,tomStatic=colorStaticTOM)
    tab<-as.data.frame(ret)
    tab$Gene=colnames(fullDatExpr)
    write.table(tab,file=paste('WGCNA_',prefix,'ClusterAssignment.tsv',sep=''),sep='\t',row.names=F)

    return(ret)
}


#' Functional enrichment
#' @param datExpr a sample by gene matrix of expression
#' @param module assignments (generally by color)
#' @return Table of most enrichmed terms for each module
getEnrichment<-function(datExpr,colorh1,prefix){
  tab<-read.table('../../data/HugoGIDsToEntrez_DAVID.txt',header=T,sep='\t',quote='"')
  eids=tab$To[match(colnames(datExpr),tab$From)]
  res=GOenrichmentAnalysis(colorh1[!is.na(eids)],eids[!is.na(eids)],organism='human',nBestP=10)
  gtab = res$bestPTerms[[4]]$enrichment
  write.table(gtab,file=paste(prefix,'top10GOTermsPermodule.csv',sep=''),sep=',')
  return(gtab)
}

#'
#'Given an expression matrix and a color assignment from the clustering
#'do some analysis
evalEigenModules<-function(datExpr,colorh1,pids=NA,prefix=''){

  datME=moduleEigengenes(datExpr,colorh1)$eigengenes
  if(is.na(pids))
    pids=rownames(datExpr)

  dissimME=(1-t(cor(datME, method="p")))/2
  hclustdatME=hclust(as.dist(dissimME), method="average" )
  # Plot the eigengene dendrogram
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based of the module eigengenes")

  ##now try to plot clusters
  which.module='orange'
  #sizeGrWindow(8,7);
 # which.module="green"
  for(which.module in unique(colorh1)){
    ME=datME[, paste("ME",which.module, sep="")]
    pdf(paste(prefix,which.module,'moduleInGeneExpression.pdf',sep=''))
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
    plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
            nrgcols=30,rlabels=F,rcols=which.module,clabels=pids,
            main=which.module, cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,
            ylab="eigengene expression",xlab="Patient Sample")
    dev.off()
}
  rownames(datME)<-rownames(datExpr)
  write.table(paste(prefix,'eigenModules.tsv',sep=''),sep='\t')
  return(datME)
}
