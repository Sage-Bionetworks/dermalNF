import os,re,sys
import pandas
import math

first_bamdir='/scratch/DERMALNF/disk1'
second_bamdir='/scratch/DERMALNF/disk3/wgs_analysis/'

first_fis=[a for a in os.listdir(first_bamdir) if '.bam' in a]
second_fis=[a for a in os.listdir(second_bamdir) if '.bam' in a]

aws_cmd='aws s3 cp'
s3_bucket='s3://ctf-data/wgs-bams/'

for file in first_fis:
    cmd=aws_cmd+' '+file+' '+s3_bucket+file
