###compare RNAseq data to CNV/proteomics data

source('../../bin/dermalNFData.R')

#mapped to gene
rna_counts=rna_count_matrix()

#mapped to gene
proteomics=prot_normalized()

#so this data isn't mapped to a gene yet
cnv=cnv_segmented()

##first let's calculate and upload matrix file
source("../../bin/clusterCNVBySample.R")
main()

sf=File(filename,parentId='syn5014748')
used(sf)<-synGet('syn5015035')
synStore(sf)