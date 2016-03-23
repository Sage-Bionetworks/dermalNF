## get TCGA expression data from synapse


library(synapseClient)
library(data.table)
library(ggplot2)
synapseLogin()

if(!exists('alldat')){
    #This is older
    #tumordata<-synGet('syn1695373')
    tumordata<-synGet('syn4311114')
    alldat<<-as.data.frame(fread(tumordata@filePath))
}
##now get mRNA patients
#source("../../bin/TcgaMutationalData.R")
#mut.data<-getNF1Data()

#mut.pats=sapply(as.character(mut.data$Tumor),function(x) paste(unlist(strsplit(x,split='-'))[1:4],collapse='-'))
#rna.pats<-sapply(colnames(alldat),function(x) paste(unlist(strsplit(x,split='-'))[1:4],collapse='-'))

#if(!exists('pat.dis'))
#  pat.dis<<-synTableQuery('SELECT acronym,samples from syn3281840')

#tumsByDis<-sapply(unique(pat.dis@values$acronym),function(x) pat.dis@values$samples[which(pat.dis@values$acronym==x)])

