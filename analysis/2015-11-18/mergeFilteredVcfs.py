'''
VCFs should be merged by sample
'''
import sys,os,re
sys.path.append('../../bin/')
import  wgsVariantRefinement as wgs

##get all synfile paths
synfiles=[wgs.syn.get(si,downloadFile=True).path+'\n' for si in wgs.filtered_syn_ids]
#now we need to gzip all files
print 'Gzipping files'
for a in synfiles:
	if not os.path.exists(a+'.gz'):
		os.system('bgzip '+a)

gzippedfiles=[re.sub('.vcf','.vcf.gz',a) for a in synfiles]
open('allVcfFiles.txt','w').writelines(gzippedfiles)
#gzippedfiles=[a.strip() for a in open('allVcfFiles.txt','r').readlines()]
#now call bcftools to merge into large vcf
vcf='allHardFilteredDermalNFSamples.vcf'

cmd='bcftools merge -o %s -l %s --use-header %s'%(vcf,'allVcfFiles.txt',gzippedfiles[0])
print cmd
#now doing bcftools merge
os.system(cmd)
#now call vqsr?
#res=wgs.snp_calibrate_model(vcf)
