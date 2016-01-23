##new plots of NF1 mutations, somatic mutations across all cancer genes

source("../../bin/WGSData.R")
nf1.res=getMutationStatsForGene('NF1')
hras.res=getMutationStatsForGene('HRAS')

##also get CNV data to process/visualize
#source("../../bin/segmentCNVData.R")

