#!/bin/bash

vcf='../2015-11-20/recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf'

genebed='../../data/UCSC_hg19_knownGene_2015_11_24.gtf'

cmd="bedtools intersect -a $vcf -b $genebed"

echo $cmd
$cmd > geneRegionsOnly.vcf

##now we have to re-add header

bcftools view -h $vcf > headerFile.txt

cat headerFile.txt geneRegionsOnly.vcf > geneRegionsWH.vcf
