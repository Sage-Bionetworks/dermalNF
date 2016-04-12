##compare proteomics to TCGA

source("../../bin/dermalNFData.R")
dermal.prots<-prot_normalized()

##now get TCGA data
#pancan.rppa<-read.table(synGet('syn1710429')@filePath,comment='',header=T,row.names=1)


pancan.rppa<-read.table('TCGA-PANCAN16-RBN.csv',sep=',',header=T)
pancan.meta<-read.table('TCGA-PANCAN16-RBN-meta.tsv',sep='\t',header=T)

gene.names<-list(Tuberin='TSC2',YB.1='YBX1',STAT5.alpha='STAT5A',
                 STAT3_pY705='STAT3',Shc_pY317='SCH1',Transglutaminase='TGM2',
                 Stathmin='STMN1',S6='RPS6',Rab.25='RAB25',PI3K.p85='PIK3R1/2',PAI.1='SERPINE1',
                 Paxillin='PXN',beta.Catenin='CTNNB1',Beclin='BECN1',Bcl2='BCL2',Claudin.7='CLDN7',
                 Collagen_VI='COL6A1',Cyclin_B1='CCNB1',Cyclin_D1='CCND1',Cyclin_E1='CCNE1',Cyclin_E2='CCNE2',
                 DJ.1='PARK7',E.Cadherin='CDH1',Rab.25='RAB25',mTOR='FRAP1',
                 MEK1_pS217_S221='MAP2K1',Heregulin='NRG1',ETS.1='ETS1',ER.alpha='ESR1',
                 Bap1.c.4='BAP1',AMPK_alpha='AMPK',Annexin.1='ANXA1',Annexin_VII='ANXA7',
                 P.Cadherin='CDH3',N.Cadherin='CDH2',Myosin.IIa_pS1943='MYH9',c.Jun_pS73='JUN',
                 c.Kit='KIT',c.Myc='MYC',C.Raf='RAF1',AMPK_alpha='PRKAA1')


matches<-sapply(dermal.prots$Protein,function(x){
  m1=grep(paste('^',x,'$',sep=''),colnames(pancan.rppa),ignore.case=T)
  if(length(m1)==0){
    m1<-grep(x,gene.names)
    if(length(m1)>0)
      return(names(gene.names)[m1])
    else
      return(NULL)
  }else{
    return(colnames(pancan.rppa)[m1])
  }
})
um<-unlist(matches)

pm<-t(pancan.rppa[,unlist(matches)])
colnames(pm)<-pancan.rppa$TCGA_patient_barcode
rownames(pm)<-names(um)
dm<-dermal.prots[match(names(um),dermal.prots$Protein),-1]
rownames(dm)<-names(um)

zvals<-which(apply(dm,2,function(x) all(x==0)))

res<-sapply(unique(colnames(pm)),function(x) if(length(which(colnames(pm)==x))>1) return(rowMeans(pm[,which(colnames(pm)==x)],na.rm=T)) else return(pm[,x]))
pc=prcomp(t(cbind(res,dm[,-zvals])))
library(ggbiplot)

tgroups<-as.character(pancan.rppa$Tumor)
names(tgroups)<-pancan.rppa$TCGA_patient_barcode
tgroups<-tgroups[unique(names(tgroups))]

dgroups<-rep("DermalNF",ncol(dm[,-zvals]))
names(dgroups)<-colnames(dm)[-zvals]

ggbiplot(pc,var.axes=FALSE,groups=c(tgroups,dgroups),ellipse=TRUE)
pc2=prcomp(t(cbind(res[,which(tgroups%in%c('LGG','GBM','COAD','READ','LUSC','LUAD','SKCM'))],dm[,-zvals])))
t2groups<-tgroups[which(tgroups%in%c('LGG','GBM','COAD','READ','LUSC','LUAD','SKCM'))]

ggbiplot(pc2,var.axes=FALSE,groups=c(t2groups,dgroups),ellipse=TRUE)




