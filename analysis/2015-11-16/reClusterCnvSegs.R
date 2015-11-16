source("../../bin/clusteRCNVBySample.R")
main()

##now put them all onto synapse
allfiles=list.files('./')
allfiles=allfiles[-grep('.R$',allfiles)]

for(a in allfiles){
  used=list(list(entity='syn5049753',wasExecuted=F),
            list(name='clusterCNVBySample.R',wasExecuted=T,url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2015-11-16/reClusterCnvSegs.R'))
  
  synStore(File(a,parentId='syn5049702'),used=used,activityName='clustering segments based on variability')

}
