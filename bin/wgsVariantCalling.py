#!/usr/local/bin/python
'''
wgsVariantCalling tries to identify differences in variants between
samples. this is done using vcfs because the BAMS have already been
processed....

https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
'''

_author_='Sara JC Gosline'
_email_='sara.gosline@sagebase.org'


import sys,os,re
import argparse

##goal is to get unfiltered vcf files from dermal nf and execute a series of
##java/r commands on each.

#get all files from wgs vcf directory
wgs_vcf='syn4984931'


import synapseclient
syn = synapseclient.Synapse()
syn.login()

query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")

#only get those that are not filtered already

syn_ids=[s['entity.id'] for s in query_res['results'] if 'hard-filtered' not in s['entity.name']]
filtered_syn_ids=[s['entity.id'] for s in query_res['results'] if 'hard-filtered' in s['entity.name']]


'''
From
 java -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R reference.fasta \
   -V hapmap.vcf \
   --discordance myCalls.vcf \
   -o output.vcf \
   -sn mySample


'''

def find_variants_between_samples(fn,libdir='../../lib/',gatkDir='../../../'):

    ref=os.path.join(libdir,'ucsc.hg19.fasta')

    base=os.path.basename(fn).split('.')[0]
    gatk=os.path.join(gatkDir,'GenomeAnalysisTK.jar')
    #ref=os.path.join(libdir,'human_g1k_v37.fasta')
    gatk_commnt='java -jar GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R reference.fasta \
   -V hapmap.vcf \
   --discordance myCalls.vcf \
   -o output.vcf \
   -sn mySample'
   os.system(gatk_command)
   return recal,tran,rscr

def snp_apply_model(fn,recalFile,trancheFile,libdir,gatkDir):
    '''
    This calls the GATK command for snps
    '''

    output='recalibrated_snps_for_'+os.path.basename(fn)
    ref=os.path.join(libdir,'ucsc.hg19.fasta')
    gatk=os.path.join(gatkDir,'GenomeAnalysisTK.jar')

    cmd='java -jar %s \
    -T ApplyRecalibration \
    -R %s \
    -input %s \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile %s \
    -tranchesFile %s \
    -o %s'%(gatk,ref,fn,recalFile,trancheFile,output)
    print(cmd)
    os.system(cmd)
    return output

'''
    Now we do the same song and dance for indels - decide on parameters,
    calibrate model apply model
    Step 5: calibrate model

    Step 6: apply model


'''

def indel_calibrate_model(fn,libdir='../../lib/',gatkDir='../../../'):
    '''
    calibrate model for indels now
    '''
    mills=os.path.join(libdir,'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf')
    ref=os.path.join(libdir,'ucsc.hg19.fasta')

    base=os.path.basename(fn).split('.')[0]
    tran='recal_INDEL_'+base+'.tranches'
    recal='recal_INDEL_'+base+'.recal'
    rscr='recalibrate_INDEL_'+base+'_plots.R'
    gatk=os.path.join(gatkDir,'GenomeAnalysisTK.jar')
    #ref=os.path.join(libdir,'human_g1k_v37.fasta')
    gatk_command='java -jar %s \
    -T VariantRecalibrator \
    -R %s \
    -input %s \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 %s \
    -an DP \
    -an QD \
    -an FS \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile %s \
    -tranchesFile %s \
    -rscriptFile %s '%(gatk,ref,fn,mills,recal,tran,rscr)
    if os.path.exists(recal) and os.path.exists(tran):
        print 'Recalibrated files already exist, will not recalculate'
        return recal,tran,rscr
    else:
    	print gatk_command

    	os.system(gatk_command)
    	return recal,tran,rscr


def indel_apply_model(fn,recalFile,trancheFile,libdir,gatkDir):
    '''
    This calls the GATK command for indels
    '''

    output='recalibrated_indels_for_'+os.path.basename(fn)
    ref=os.path.join(libdir,'ucsc.hg19.fasta')
    gatk=os.path.join(gatkDir,'GenomeAnalysisTK.jar')

    cmd='java -jar %s \
    -T ApplyRecalibration \
    -R %s \
    -input %s \
    -mode INDEL \
    --ts_filter_level 99.0 \
    -recalFile %s \
    -tranchesFile %s \
    -o %s'%(gatk,ref,fn,recalFile,trancheFile,output)
    print(cmd)
    os.system(cmd)
    return output


def main():
    parser=argparse.ArgumentParser(description='Run GATK using the command line')
    ##we can take either a vcf or synapse id
    parser.add_argument('--vcf',dest='vcf_file',help='Path to VCF file to process')
    parser.add_argument('--synapseId',dest='synapse_id',help='Synapse ID of VCF to process')
    ##paths to library and GATK
    parser.add_argument('--libDir',dest='libdir',help='Relative path the reference fasta and vcf files')
    parser.add_argument('--gatkDir',dest='gatkdir',help='Relative path to GATK')
    #model type
    parser.add_argument('--modelType',dest='model',default='snp',help='Either SNP,INDEL both, depending on what type of model to calibrate/apply')

    args=parser.parse_args()

    if args.vcf_file is None:
        if args.synapse_id is None:
            print "Need to have either a VCF or synapse identifier, try again"
            sys.exit()
        else:
            vcf=synGet(args.synapse_id).path
    else:
        vcf=args.vcf_file

    if args.libdir is None or args.gatkdir is None:
        print 'Need to have both library directory and GATK directory to run!'
        sys.exit()

    if args.model.lower()=='snp' or args.model.lower()=='both':
        print 'Calibrating and applying SNP model to %s'%(vcf)
        recal,tran,rscr = snp_calibrate_model(vcf,args.libdir,args.gatkdir)
        print 'Finished recalibrating model, run %s for more analysis.  Now  applying model'%(rscr)
        outputfile=snp_apply_model(vcf,recal,tran,args.libdir,args.gatkdir)
        print 'Finished applying model, see'+outputfile
        vcf=outputfile

    if args.model.lower()=='indel' or args.model.lower()=='both':
        print 'Calibrating and applying INDEL model to %s'%(vcf)
        recal,tran,rscr = indel_calibrate_model(vcf,args.libdir,args.gatkdir)
        print 'Finished recalibrating model, run %s for more analysis.  Now applying model'%(rscr)
        outputfile=indel_apply_model(vcf,recal,tran,args.libdir,args.gatkdir)
        print 'Finished applying model, see'+outputfile



            #now for each synapse id, get file and runcode on it!
#    for si in syn_ids:
#        so=syn.get(si)
#        vcf=so.path
if __name__=='__main__':
     main()
