'''
Now that we have filtered VCF we can create MAF for each tumor vs. normal
'''

vcf2maf='perl ../../../vcf2maf-master/vcf2maf.pl'

##test out one tumor id and one normal id for a single patient
wgs_vcf='syn4984931'


import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()
#get all files for patient1
query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")

##now get metadata for all files
syn_files=dict([(res['entity.id'],res['entity.name']) for res in query_res['results'] if 'hard-filtered' not in res['entity.name']])
syn_ids=syn_files.keys()

##now get metadata only for these files to determine which synapse ids belong to
##a specific patient
annotes=[syn.getAnnotations(a) for a in syn_ids]
reffasta='../../lib/ucsc.hg19.fasta'

for pat in ['1','2','3','4','5','6','7','8','9','10','11','12','13']:

    p1ids=[a for a in annotes if a['patientID'][0] in ['CT00000'+pat,'CT0000'+pat]]

    blood=[a['id'] for a in p1ids if a['tissueID'][0]=='PBMC']
    tumor=[a['id'] for a in p1ids if a['tissueID'][0]!='PBMC']
    print 'Found %d tumor samples and %d blood samples for patient %s'%(len(tumor),len(blood),pat)

    vcffile='recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf.gz'

    if len(blood)>0 and blood[0] in syn_files.keys():
        bloodfile=re.sub('.vcf','',syn_files[blood[0]])

    cmfile='patient_'+pat+'_commands.sh'
    patsh=open(cmfile,'w')
    for t in tumor:
        tumfile=re.sub('.vcf','',syn_files[t])
        print t,blood[0]
        ##first create vcf file with only two samples
        bcftoolscmd=''
        #then create a symlink from that file to a 'ghost' .vcf file

        #then run vcf2maf on that subset
        cmd=vcf2maf+" --input-vcf ../2015-12-02/recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf --tumor-id %s --normal-id %s --output-maf tumorVsNormal_pat%s_%s.maf \
        --vcf-tumor-id %s --vep-forks 8 --species homo_sapiens \
        --vcf-normal-id %s --ref-fasta %s"%(tumfile,bloodfile,pat,t,tumfile,bloodfile,reffasta)
        print cmd
        patsh.write(cmd+'\n')
        #os.system(cmd)
    patsh.close()
    os.system('sh '+cmfile+' &')

        #examples for patient 1
'''
syn4987466 syn4985519
syn4985596 syn4985519
syn4985378 syn4985519
##here is the test command, as far as I can tell...
perl ../../../vcf2maf-master/vcf2maf.pl --input-vcf recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf --tumor-id=syn4987466 --normal-id=syn4985519 --output-maf=tumorVsNormal_pat1.maf
'''
