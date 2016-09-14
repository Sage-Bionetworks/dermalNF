require(synapseClient)
synapseLogin()

orig.dir<-'/scratch/DERMALNF/disk1/proteomics_data'

for(file in list.files(orig.dir)){
    tab<-read.table(paste(orig.dir,file,sep='/'),sep='\t',header=T)
    orig<-tab[,1:5]
    for(i in colnames(tab)[6:9]){
        nzvals<-which(tab[,1]!=0.0)
        print(paste('found',length(nzvals),'non zero values for column',i))
        newtab<-cbind(orig[nzvals,],tab[nzvals,i])
        fname=gsub('.txt',paste(gsub('.','_',i,fixed=T),'.txt',sep=''),file)
        write.table(newtab,file=fname)
    }
}
