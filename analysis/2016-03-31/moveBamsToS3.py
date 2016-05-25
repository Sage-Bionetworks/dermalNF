import os,re,sys

first_bamdir='/scratch/DERMALNF/disk1/'
second_bamdir='/scratch/DERMALNF/disk3/wgs_analysis/'

first_fis=[a for a in os.listdir(first_bamdir) if '.bam' in a]
second_fis=[a for a in os.listdir(second_bamdir) if '.bam' in a]

aws_cmd='aws s3 cp'
s3_bucket='s3://ctf-data/wgs-bams/'

fpath=first_bamdir#+file
cmd=aws_cmd+' '+fpath+' '+s3_bucket+' --recursive --exclude "*" --include "*.bam"'
print cmd
os.system(cmd+' &')

#for file in second_fis:
fpath=second_bamdir#+file
cmd=aws_cmd+' '+fpath+' '+s3_bucket+' --recursive --exclude "*" --include "*.bam"'
print cmd
os.system(cmd+' &')
