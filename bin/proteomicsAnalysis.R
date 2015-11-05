##use this file to process the iTRAQ data
require(synapseClient)
synapseLogin()
#nfdat<-synGet('syn4984627')


all.prots<-NULL
for(i in 1:ncol(res))
  all.prots<-union(all.prots,res[['Prot.ids',i]])

expr.prots<-res[['Prot.ids',1]]
for(i in 2:ncol(res))
  expr.prots<-intersect(expr.prots,res[['Prot.ids',i]])


prot.ids<-unique(unlist(sapply(all.prots,function(x) unlist(strsplit(x,split=';')))))

#now create biomart mapping
require(biomaRt)
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

epep="ensembl_peptide_id"
egene='hgnc_symbol'

gene.mapping<-getBM(attributes=c(epep,egene),filters=c(epep),values=as.list(prot.ids),mart=ensembl)

#now we have ratios and raw values, and controls - let's create heatmap of each

ratio.mat<-sapply(all.prots,function(x){
  pvec<-NULL
  samps<-NULL
  for(i in 1:ncol(res)){
      rv=grep(x,res[['Prot.ids',i]])
      if(length(rv)==0)
          pvec<-c(pvec,rep(0,3))
      else
          pvec<-c(pvec,res[[1,i]][rv,])
    samps<-c(samps,colnames(res[[1,i]]))
  }
  names(pvec)<-samps
  unlist(pvec)
})



expr.ratio.mat<-sapply(expr.prots,function(x){
  pvec<-NULL
  samps<-NULL
  for(i in 1:ncol(res)){
      rv=grep(x,res[['Prot.ids',i]])
      if(length(rv)==0)
          pvec<-c(pvec,rep(0,3))
      else
          pvec<-c(pvec,res[[1,i]][rv,])
    samps<-c(samps,colnames(res[[1,i]]))
  }
  names(pvec)<-samps
  unlist(pvec)
})

expr.ratio.mat[which(is.na(expr.ratio.mat),arr.ind=T)]<-0.0
expr.ratio.mat<-expr.ratio.mat[-grep('EMPTY',rownames(expr.ratio.mat)),]
require(pheatmap)

#now remove zeroes
expr.ratio.mat<-expr.ratio.mat[,-which(apply(expr.ratio.mat,2,function(x) mean(x,na.rm=T))==0)]

samp.mapping<-sapply(rownames(expr.ratio.mat),function(x){
    for(i in 1:ncol(res))
        if(x%in%colnames(res[[1,i]]))
            return(colnames(res)[i])
})


new.cn<-sapply(colnames(expr.ratio.mat),function(x) if(x%in%gene.mapping[,1]) return(gene.mapping[which(gene.mapping[,1]==x),2]) else return(x))
gene.expr.ratio.mat<-expr.ratio.mat
colnames(gene.expr.ratio.mat)<-new.cn
pheatmap(log2(as.matrix(gene.expr.ratio.mat)+0.001),cellheight=10,cellwidth=10,filename='dermalNF_topProts_allsamps_genenames.pdf',annotation_row=data.frame(Sample=samp.mapping))
############################################

                                        #plot only those non-zero values
red.mat<-gene.expr.ratio.mat[,-which(apply(gene.expr.ratio.mat,2,function(x) any(x==0)))]
pheatmap(log2(red.mat+0.001),cellheight=10,cellwidth=10,filename='dermalNF_topProts_allsamps_nonzero_genenames.pdf',annotation_row=data.frame(Sample=samp.mapping))

#############################################

                                        #now we can plot raw values
expr.raw.mat<-sapply(expr.prots,function(x){
  pvec<-NULL
  samps<-NULL
  for(i in 1:ncol(res)){
      rv=grep(x,res[['Prot.ids',i]])
      if(length(rv)==0)
          pvec<-c(pvec,rep(0,3))
      else
          pvec<-c(pvec,res[[2,i]][rv,])
    samps<-c(samps,colnames(res[[2,i]]))
  }
  names(pvec)<-samps
  unlist(pvec)
})


expr.raw.mat[which(is.na(expr.raw.mat),arr.ind=T)]<-0.0
expr.raw.mat<-expr.raw.mat[-grep('EMPTY',rownames(expr.raw.mat)),]
require(pheatmap)

#now remove zeroes
expr.raw.mat<-expr.raw.mat[,-which(apply(expr.raw.mat,2,function(x) mean(x,na.rm=T))==0)]

samp.mapping<-sapply(rownames(expr.raw.mat),function(x){
    for(i in 1:ncol(res))
        if(x%in%colnames(res[[1,i]]))
            return(colnames(res)[i])
})

##also can we cluster just controls together?

new.cn<-sapply(colnames(expr.raw.mat),function(x) if(x%in%gene.mapping[,1]) return(gene.mapping[which(gene.mapping[,1]==x),2]) else return(x))
gene.expr.raw.mat<-expr.raw.mat
colnames(gene.expr.raw.mat)<-new.cn
pheatmap(log10(as.matrix(gene.expr.raw.mat)+0.001),cellheight=10,cellwidth=10,filename='dermalNF_topProts_raw_allsamps_genenames.pdf',annotation_row=data.frame(Sample=samp.mapping))

#############################################

                                        #now do control values
expr.ctrl.mat<-sapply(expr.prots,function(x){
  pvec<-NULL
  samps<-NULL
  for(i in 1:ncol(res)){
      rv=grep(x,res[['Prot.ids',i]])
      if(length(rv)==0)
          pvec<-c(pvec,0)
      else
          pvec<-c(pvec,res[[3,i]][rv])
    samps<-c(samps,colnames(res)[i])
  }
  names(pvec)<-samps
  unlist(pvec)
})


expr.ctrl.mat[which(is.na(expr.ctrl.mat),arr.ind=T)]<-0.0

require(pheatmap)

#now remove zeroes
expr.ctrl.mat<-expr.ctrl.mat[,-which(apply(expr.ctrl.mat,2,function(x) mean(x,na.rm=T))==0)]

#pdf('dermalNF_topProts_allsamps.pdf')
pheatmap(as.matrix(expr.ctrl.mat),cellheight=10,cellwidth=10,filename='dermalNF_topProts_allsamps.pdf')
#dev.off()

##also can we cluster just controls together?

new.cn<-sapply(colnames(expr.ctrl.mat),function(x) if(x%in%gene.mapping[,1]) return(gene.mapping[which(gene.mapping[,1]==x),2]) else return(x))
gene.expr.ctrl.mat<-expr.ctrl.mat
colnames(gene.expr.ctrl.mat)<-new.cn
pheatmap(log10(as.matrix(gene.expr.ctrl.mat)+0.001),cellheight=10,cellwidth=10,filename='dermalNF_topProts_ctrl_vals_genenames.pdf')
