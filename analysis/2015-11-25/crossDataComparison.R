###compare RNAseq data to CNV/proteomics data

source('../../bin/dermalNFData.R')

#mapped to gene
rna_counts=rna_count_matrix()

#mapped to gene
proteomics=prot_normalized()

#so this data isn't mapped to a gene yet
#cnv=cnv_segmented()

##first let's calculate and upload matrix file
source("../../bin/clusterCNVBySample.R")
main()

filename='medianlogRRatioValuesBygene.txt'
sf=File(filename,parentId='syn5049702')
#act=Activity(name='built matrix of median segments by gene and region',description='main() method from clusterCNVBySample.R')
used(sf)<-list(list(name='crossDataComparison.R',url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2015-11-25/crossDataComparison.R'),synGet('syn5049753'))
synStore(sf,activityName='built matrix of median segments by gene')

filename='medianlogRRatioValuesByregion.txt'
sf=File(filename,parentId='syn5049702')
#act=Activity(name='built matrix of median segments by gene and region',description='main() method from clusterCNVBySample.R')
used(sf)<-list(list(name='crossDataComparison.R',url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2015-11-25/crossDataComparison.R'),synGet('syn5049753'))
synStore(sf,activityName='built matrix of median segments by region')

##now we can use the snp data per gene using a single function call
source('../../bin/dermalNFData.R')

cnv<-cnv_segmented_by_gene()

##last we need to get the mappings from one gene to another. this will be in the annotations...
rna.pat.tum.ids<-list()
other.pat.tum.ids<-list()

#now we can try to map one sample to the other...
