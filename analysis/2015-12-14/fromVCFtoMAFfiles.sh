#!/bin/bash

merged_vcf=''
gdir='../../../'
libdir='../../lib/'
cmd='python ../../bin/wgsVariantRefinement.py --vcf $merged_vcf --libDir $libdir --gatkDir $gdir --modelType both'
echo $cmd
$cmd

newvcf=''

