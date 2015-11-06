#!/usr/local/bin/python
'''
wgsVariantRefinement follos the GATK toolkit starting with our original vcfs
coming from CTF. We skip a few steps and go straight to variant refinement here:

https://www.broadinstitute.org/gatk/guide/article?id=2805
'''

_author_='Sara JC Gosline'
_email_='sara.gosline@sagebase.org'


import sys,os,re


##goal is to get unfiltered vcf files from dermal nf and execute a series of
##java/r commands on each.

#get all files from wgs vcf directory
wgs_vcf='syn4984931'


import synapseclient
syn = synapseclient.Synapse()

query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")

#only get those that are not filtered already
syn_ids=[s['entity.id'] for s in query_res['results'] if 'hard-filtered' not in s['entity.name']]

'''
From

Step 2: build SNP recalibration model
java -jar GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R reference.fa \
    -input raw_variants.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an InbreedingCoeff \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R

Step 3: apply recalibration to data
java -jar GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R reference.fa \
    -input raw_variants.vcf \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o recalibrated_snps_raw_indels.vcf

'''

def snp_calibrate_model(fn):
    gatk_command='java -jar GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R reference.fa \
    -input %s \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -an InbreedingCoeff \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R '%(fn)
    os.sys(gatk_command)


def snp_apply_model(fn):
    '''
    This calls the GATK command for indels
    '''

'''
    Now we do the same song and dance for indels - decide on parameters,
    calibrate model apply model
    Step 5: calibrate model

    Step 6: apply model


'''

def indel_calibrate_model(fn):
    '''
    calibrate model for indels
    '''

def indel_apply_model(fn):
    '''
    apply indel model to data
    '''


#now for each synapse id, get file and runcode on it!
for si in syn_ids:
    so=syn.get(si)
    vcf=so.path
