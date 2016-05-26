source("../../dermalNF/bin/WGSData_VarDict.R")


for(syn in c('syn6047991','syn6087182')){
  
  all.gene.muts<-read.table(synGet(syn)@filePath,sep='\t',header=T)


  transcripts<-read.table(synGet('syn5579597')@filePath,header=T)

  #from mutant data, create mappin gof all genes and transcripts identified....
  gene.mut<-data.frame(Transcript=unique(as.character(all.gene.muts$Transcript)),
                     shortened=sapply(unique(as.character(all.gene.muts$Transcript)),function(x) unlist(strsplit(x,split='.',fixed=T))[1]))
  gene.mut$Gene=all.gene.muts$Gene[match(gene.mut$Transcript,all.gene.muts$Transcript)]


  
  #now do the same for all RNA seq data
  gene.expr<-data.frame(Gene=sapply(rownames(transcripts),function(x) unlist(strsplit(x,split='.',fixed=T))[1]),
                        Transcript=sapply(rownames(transcripts),function(x) unlist(strsplit(x,split='.',fixed=T))[2]),
                        shortened=sapply(rownames(transcripts),function(x) paste(unlist(strsplit(unlist(strsplit(x,split='.',fixed=T))[2],split='_'))[1:2],collapse='_')))
  
  gene.expr$Present<-apply(transcripts,1,function(x) any(x>1))
  
  
  overlap<-intersect(gene.expr$shortened,gene.mut$shortened)
  
  #here are the genes covered, not all are unique! :()
  genes.covered=unique(data.frame(Mut=gene.mut$Gene[match(overlap,gene.mut$shortened)],Exp=gene.expr$Gene[match(overlap,gene.expr$shortened)]))
  
  
  ##genes for which there are mutations but no expression
  missing=setdiff(gene.mut$shortened,gene.expr$shortened)
  
  ##are they covered by other gene? 
  genes.missed<-unique(gene.mut$Gene[match(missing,gene.mut$shortened)])
  
  ##here are the genes are ain the mutation data and not covered in the expression data...
  not.covered=setdiff(genes.missed,gene.expr$Gene)
  
  mut.inds<-match(overlap,gene.mut$shortened)
  expr.inds<-match(overlap,gene.expr$shortened)
  
  combined<-cbind(Expr=gene.expr[expr.inds,],Mut=gene.mut[mut.inds,])
  
  expr.muts<-combined$Mut.Transcript[which(combined$Expr.Present)]
  
  expr.gene.muts<-subset(all.gene.muts,Transcript%in%expr.muts)
  
  fname=gsub('.tsv','_filteredForExpr.tsv',basename(synGet(syn)@filePath))
  write.table(expr.gene.muts,file=fname,sep='\t',quote=F)
  
  synStore(File(fname,parentId='syn5605256'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-17/filterByTranscripts.R')),used=list(list(entity='syn6047991')))
}

