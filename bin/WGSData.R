##process WGS data

library(synapseClient)
synapseLogin()

allmafs<-synapseQuery("select * from entity where parentId=='syn5522808'")

getMutationSummary<-function(){

    allMuts<-sapply(allmafs$entity.id,function(x){
        res<-synGet(x)
        tab<-read.table(gzfile(res@filePath),sep='\t',header=T)

    })
}
#summary(tab$Consequence)

library(parallel)

getMutationMatrix<-function(impact='HIGH'){

    allMuts<-mclapply(allmafs$entity.id,function(x){
        res<-synGet(x)
        tab<-read.table(gzfile(res@filePath),sep='\t',header=T)
        #dont filter by consequence
        #vars<-tab[which(tab$Consequence%in%mutClasses),]
        if(!is.na(impact))
            vars<-tab[which(tab$IMPACT==impact),]
        else
            vars<-tab
        print(length(which(vars$Hugo_Symbol)=='NF1'))#summary(vars$Hugo_Symbol))
        return (vars)
    },mc.cores=4)
}
