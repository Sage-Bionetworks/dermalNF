
source("../../bin/WGSData.R")

synStore(File('v1.3_anno_191025841cf776d395abb257495fc835.tsv',parentId='syn5605256',name='CADDScoresForNF1Region'))

tab<-read.table(synGet('syn5820042')@filePath,comment.char='',skip=1,header=T)

rtab<-tab[-c(which(is.na(tab$GeneName)),grep("CTD",tab$GeneName),grep("RP1",tab$GeneName),grep("RN",tab$GeneName),grep('CTC',tab$GeneName),grep('AC',tab$GeneName)),]
require(ggplot2)

p<-ggplot(subset(rtab,PHRED>10))+geom_histogram(aes(x=PHRED,fill=GeneName),position='stack')
print(p+ggtitle("CADD scores over 10 in chr17 region"))




pn<-ggplot(subset(rtab,GeneName=='NF1'))+geom_point(aes(x=Pos,y=PHRED,colour=Consequence,size=2))

pa<-ggplot(subset(rtab,PHRED>15))+geom_point(aes(x=Pos,y=PHRED,colour=Consequence,size=2))

nf1.deets=subset(rtab,GeneName=='NF1')[,c('Pos','PHRED','Consequence')]

nf1.muts<-getMutationStatsForGene('NF1')

##first get idx values
mut.idx<-match(nf1.muts$Start_Position,nf1.deets$Pos)
na.vals<-which(is.na(mut.idx))
mut.idx[na.vals]<-sapply(na.vals,function(x){
  match(nf1.muts[x,'Start_Position']-1,nf1.deets$Pos)
})

nf1.muts$PHRED<-nf1.deets$PHRED[mut.idx]
nf1.muts$Consequence<-nf1.deets$Consequence[mut.idx]
mut.score=nf1.muts$PHRED
names(mut.score)<-nf1.muts$Protein_Change

files<-heatmapWithPhredScore(nf1.muts,'nf1Mutation.png',nf1.deets,10)
for(f in files){
  sf=File(f,parentId='syn5605256')
  this.script='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-22/processCaddScores.R'
  synStore(sf,used=list(list(entity='syn5820042'),list(url=this.script,wasExecuted=TRUE)))
  
}

