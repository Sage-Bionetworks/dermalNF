##plot germline/somatic mutations from snpEff/varDict

tab<-read.table('nf1_combined.maf',header=T,as.is=T,quote='"',sep='\t')

passed=subset(tab,PASS=="TRUE")

res=unique(tab[which(tab$Effect%in%c('MODERATE','HIGH')),c('Sample','Start','Amino_Acid_Change','cDNA_Change','Status','Type')])
View(res)