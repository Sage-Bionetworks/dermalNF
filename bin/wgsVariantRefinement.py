#!/usr/local/bin/python
'''
wgsVariantRefinement follos the GATK toolkit starting with our original vcfs
coming from CTF. We skip a few steps and go straight to variant refinement here:

https://www.broadinstitute.org/gatk/guide/article?id=2805
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

def snp_calibrate_model(fn,libdir='../../lib/',gatkDir='../../../'):
    hapmap=os.path.join(libdir,'hapmap_3.3.hg19.sites.vcf')
    omni=os.path.join(libdir,'1000G_omni2.5.hg19.sites.vcf')
    otg=os.path.join(libdir,'1000G_phase1.snps.high_confidence.hg19.sites.vcf')
    dbsnp=os.path.join(libdir,'dbsnp_138.hg19.vcf')
    ref=os.path.join(libdir,'ucsc.hg19.fasta')

    base=os.path.basename(fn).split('.')[0]
    tran='recal_SNP_'+base+'.tranches'
    recal='recal_SNP_'+base+'.recal'
    rscr='recalibrate_SNP_'+base+'_plots.R'
    gatk=os.path.join(gatkDir,'GenomeAnalysisTK.jar')
    #ref=os.path.join(libdir,'human_g1k_v37.fasta')
    gatk_command='java -jar %s \
    -T VariantRecalibrator \
    -R %s \
    -input %s \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 %s  \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 %s \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s \
    -an DP \
    -an QD \
    -an FS \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile %s \
    -tranchesFile %s \
    -rscriptFile %s '%(gatk,ref,fn,hapmap,omni,otg,dbsnp,recal,tran,rscr)
    if os.path.exists(recal) and os.path.exists(tran):
        print 'Recalibrated files already exist, will not recalculate'
        return recal,tran,rscr
    else:
    	print gatk_command

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
    parser.add_argument('--modelType',dest='model',default='snp',help='Either SNP or INDEL, depending on what type of model to calibrate/apply')

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

    if args.model.lower()=='snp':
        print 'Calibrating and applying SNP model to %s'%(vcf)
        recal,tran,rscr = snp_calibrate_model(vcf,args.libdir,args.gatkdir)
        print 'Finished recalibrating model, run %s for more analysis.  Now  applying model'%(rscr)
        outputfile=snp_apply_model(vcf,recal,tran,args.libdir,args.gatkdir)
        print 'Finished applying model, see'+outputfile

    elif  args.model.lower()=='indel':
        print 'Calibrating and applying INDEL model to %s'%(vcf)
        recal,tran,rscr = indel_calibrate_model(vcf,args.libdir,args.gatkdir)
        print 'Finished recalibrating model, run %s for more analysis.  Now applying model'%(rscr)
        outputfile=indel_apply_model(vcf,recal,tran,args.libdir,args.gatkdir)
        print 'Finished applying model, see'+outputfile



            #now for each synapse id, get file and runcode on it!
#    for si in syn_ids:
#        so=syn.get(si)
#        vcf=so.path
main()
