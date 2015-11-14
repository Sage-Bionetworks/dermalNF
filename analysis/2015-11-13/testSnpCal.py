import sys,os,re
sys.path.append('../../bin/')
import  wgsVariantRefinement as wgs

vcf='allDermalNFSamples.vcf'
res=wgs.snp_calibrate_model(vcf)

