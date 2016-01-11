'''
Now that we have filtered VCF we can create MAF for each tumor vs. normal
'''

vcf2maf='perl ../../../vcf2maf-master/vcf2maf.pl'

##test out one tumor id and one normal id for a single patient
wgs_vcf='syn4984931'


import synapseclient
syn = synapseclient.Synapse()
syn.login()
#get all files for patient1
query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")

##now get metadata for all files
syn_ids=[res['entity.id'] for res in query_res['results'] if 'hard-filtered' not in res['entity.name']]

##now get metadata only for these files to determine which synapse ids belong to
##a specific patient
annotes=[syn.getAnnotations(a) for a in syn_ids]


for pat in ['1','2','3','4','5','6','7','8','9','10','11',]:

    p1ids=[a for a in annotes if a['patientID'][0] in ['CT00000'+pat,'CT0000'+pat]]

    blood=[a['id'] for a in p1ids if a['tissueID'][0]=='PBMC']
    tumor=[a['id'] for a in p1ids if a['tissueID'][0]!='PBMC']
    print 'Found %d tumor samples and %d blood samples for patient %s'%(len(tumor),len(blood),pat)

    vcffile='recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf.gz'


    for t in tumor:
        print t,blood[0]
        cmd=vcf2maf+" --input-vcf recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf \
        --tumor-id=%s --normal-id=%s --output-maf=tumorVsNormal_pat%s_%s.maf \
        --vcf-tumor-id=%s \
        --vcf-normal-id=%s"%s(t,blood[0],pat,t,t.upper(),blood[0].upper)
        print cmd
      #  os.system(cmd)

        #examples for patient 1
'''
syn4987466 syn4985519
syn4985596 syn4985519
syn4985378 syn4985519
##here is the test command, as far as I can tell...
perl ../../../vcf2maf-master/vcf2maf.pl --input-vcf recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf --tumor-id=syn4987466 --normal-id=syn4985519 --output-maf=tumorVsNormal_pat1.maf
'''
