##new plots of NF1 mutations, somatic mutations across all cancer genes

source("../../bin/WGSData.R")
nf1.res=getMutationStatsForGene('NF1')
hras.res=getMutationStatsForGene('HRAS')
tp53.res=getMutationStatsForGene('TP53')
braf.res=getMutationStatsForGene('BRAF')


allfiles=list.files('./')
ufiles=allfiles[c(grep('.png',allfiles),grep(".pdf",allfiles),grep('tsv',allfiles))]
sapply(ufiles,function(x) {
  f=File(x,parentId='syn5605256')
  synStore(f,executed=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/WGSData.R'),list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-01-23/nf1MutAnalysis.R')))
                
})

##also get CNV data to process/visualize
#source("../../bin/segmentCNVData.R")

