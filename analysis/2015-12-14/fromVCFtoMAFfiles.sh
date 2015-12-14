#!/bin/bash

merged_vcf=''
gdir='../../../'
libdir='../../lib/'
cmd='python ../../bin/wgsVariantRefinement.py --vcf $merged_vcf --libDir $libdir --gatkDir $gdir --modelType both'
echo $cmd
$cmd

zc="bgzip recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf"
echo $zc
$zc

newvcf=recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf.gz

syncmd="synapse store $newvcf --parentId syn5522790"
echo $syncmd
$syncmd

mafcmd="python ../../bin/vcfDataProcessing.py"
echo $mafcmd
$mafcmd
