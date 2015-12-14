import subprocess
import os
import argparse
import vcf


def sortGATKVCF(vcfFile, faiFile='/home/tyu/reference/GATK/ucsc.hg19.fasta.fai',sortingFunctionDir='/home/tyu/software/dermalNF/bin'):
	tempFile = vcfFile.replace(".vcf","_temp.vcf")
	sortedvcf = vcfFile.replace(".vcf","_sort.vcf")
	cmds = ['perl',os.path.join(sortingFunctionDir,'sortGATKVCF.pl'),vcfFile,faiFile,'>',tempFile]
	p = subprocess.Popen(cmds,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	stdout,stderr = p.communicate()
	if not stderr:
		print 'Sorted VCF %s -> %s' % (os.path.basename(vcfFile),os.path.basename(sortedvcf))
	else:
		print stderr
	vcf_reader = vcf.Reader(filename=vcfFile)
	vcf_writer = vcf.Writer(open(sortedvcf, 'w'), vcf_reader)
	for record in vcf_reader:
		if record.REF not in record.ALT:
			vcf_writer.write_record(record)
	os.remove(tempFile)
	return(sortedvcf)


parser = argparse.ArgumentParser()


parser.add_argument('-i', '--INPUT',dest='vcf',type=str,required = True,help='Input VCF')
parser.add_argument('-fai', '--faiFile',  dest = 'fai',type=str,required = False,help='Fai file from GATK')

args = parser.parse_args()
sortGATKVCF(args.vcf)
