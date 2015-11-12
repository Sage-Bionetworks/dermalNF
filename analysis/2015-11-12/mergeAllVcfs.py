'''
VCFs should be merged by sample
'''
import sys,os,re
sys.path.append('../../bin/')
import  wgsVariantRefinement as wgs

##get all synfile paths
synfiles=[wgs.syn.get(si).path+'\n' for si in wgs.syn_ids]
open('allVcfFiles.txt','w').writelines(synfiles)

#now call bcftools to merge into large vcf
vcf='allDermalNFSamples.vcf'

cmd='bcftools merge -o %s -l %s --use-header %s'%(vcf,'allVcfFiles.txt',synfiles[0])
print cmd

os.system(cmd)
#now call vqsr?
#res=wgs.snp_calibrate_model(vcf)
