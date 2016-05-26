#!/bin/bash

dictfile='ucsc.hg19.dict'
dbsnp='dbsnp_138.hg19.vcf'
cosmic='b37_cosmic_v54_120711_withchr.vcf'

libdir='../../lib'

cmd1="perl vcfsorter.pl $libdir/$dictfile $libdir/$dbsnp" #> $libdir/sorted_$dbsnp 2>STDERR"
echo $cmd1
$cmd1>"$libdir/sorted_$dbsnp"

cmd1="perl vcfsorter.pl $libdir/$dictfile $libdir/$cosmic" #> $libdir/sorted_$cosmic 2>STDERR"
echo $cmd1
$cmd1>"$libdir/sorted_$cosmic"
