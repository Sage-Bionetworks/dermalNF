#do some analysis for the landscape paper
source("../../bin/dermalNFData.R")

require(dplyr)
require(ggplot2)

#first do somatic load 
all.mutations<-read.table(synGet('syn5713423')@filePath,header=T,sep='\t')
cancer.gene.muts=read.table(synGet('syn5611520')@filePath,header=T,sep='\t')
non.silent=subset(all.mutations,Mutation_Type!='Silent')

som.counts=subset(non.silent,Mutation_Status=='Somatic') %>% group_by(Sample_ID) %>% summarize(MutatedGenes=n_distinct(Hugo_Symbol),DistinctMutations=n_distinct(Protein_Change))
df=tidyr::gather(som.counts,'SomaticEvent','Count',2:3)
p<-ggplot(df)+geom_bar(aes(x=Sample_ID,y=Count,fill=SomaticEvent),stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
som.cancer=subset(subset(non.silent,Hugo_Symbol%in%cancer.gene.muts$Hugo_Symbol),Mutation_Status=='Somatic') %>% group_by(Sample_ID) %>% summarize(MutatedCancerGenes=n_distinct(Hugo_Symbol),DistinctCancerMutations=n_distinct(Protein_Change))
df2=tidyr::gather(som.cancer,'SomaticEvent','Count',2:3)

##now do the same for cancer.genes?
p<-ggplot(rbind(df,df2))+geom_bar(aes(x=Sample_ID,y=Count,fill=SomaticEvent),stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+scale_y_log10()
png('somaticMutationLoadAcrossSamples.png')
print(p)
dev.off()

synStore(File('somaticMutationLoadAcrossSamples.png',parentId='syn5605256'),used=list(list(url=''),wasExecuted=TRUE))

