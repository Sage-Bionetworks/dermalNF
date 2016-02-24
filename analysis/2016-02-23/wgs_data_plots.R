##wgs plot experiments
source("../../bin/WGSData.R")

nf1.muts=getMutationStatsForGene(gene='NF1')

##here are examples of individual genes
heatmapFromMutDf(nf1.muts,fname='NF1mutations.png')

heatmapFromMutDf(getMutationStatsForGene(gene='TP53'),fname='TP53mutations.png')

heatmapFromMutDf(getMutationStatsForGene(gene='MNX1'),fname='MNX1mutations.png')


##now we want to formalize the germline mutational barrier for cancer genes