##wgs plot experiments
source("../../bin/WGSData.R")

nf1.muts=getMutationStatsForGene(gene='NF1')

##here are examples of individual genes
heatmapFromMutDf(nf1.muts,fname='NF1mutations.png')

heatmapFromMutDf(getMutationStatsForGene(gene='TP53'),fname='TP53mutations.png')

heatmapFromMutDf(getMutationStatsForGene(gene='MNX1'),fname='MNX1mutations.png')


##now we want to formalize the germline mutational barrier for cancer genes
panPatientPlots(mutTable=cancer.gene.muts,minSamples=2,notIncluded=c('Silent'))
panPatientPlots(mutTable=cancer.gene.muts,minSamples=5,notIncluded=c('Silent'))
panPatientPlots(mutTable=cancer.gene.muts,minSamples=8,notIncluded=c('Silent'))

heatmapFromMutDf(getMutationStatsForGene(gene='NUTM2A'),fname='NUTM2Amutations.png')
heatmapFromMutDf(getMutationStatsForGene(gene='KMT2C'),fname='KMT2Cmutations.png')
heatmapFromMutDf(getMutationStatsForGene(gene='NOTCH2'),fname='NOTCH2mutations.png')
