'''
This is the primary script that can be used to process the VCFs after GTK is run
on them. It requires an installation of vcf2maf.pl and VEP.
'''

##VCF2MAF installation:
vcf2maf='perl ../../../vcf2maf-master/vcf2maf.pl'
##REFERENCE FASTA installation
reffasta='../../lib/ucsc.hg19.fasta'
##OUTPUT FROM GATK:



import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()

#vcf_file=syn.get('syn5555584').path


##get all the VCF annotations so that we can process the merged file
#wgs_vcf='syn5522788'
#get all files
#query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")
##now get metadata for all files
#syn_files=dict([(res['entity.id'],res['entity.name']) for res in query_res['results'] if 'hard-filtered' not in res['entity.name']])
#syn_ids=syn_files.keys()

##now get metadata only for these files to determine which synapse ids belong to
##a specific patient
#annotes=[syn.getAnnotations(a) for a in syn_ids]

#for pat in ['1','2','3','4','5','6','7','8','9','10','11','12','13']:

    #now for each patient, determine which samples are in blood vs. tumor
 #   p1ids=[a for a in annotes if a['patientID'][0] in ['CT00000'+pat,'CT0000'+pat]]
 #   blood=[a['id'] for a in p1ids if a['tissueID'][0]=='PBMC']
 #   tumor=[a['id'] for a in p1ids if a['tissueID'][0]!='PBMC']

 #   print 'Found %d tumor samples and %d blood samples for patient %s'%(len(tumor),len(blood),pat)

 #   if len(blood)>0 and blood[0] in syn_files.keys():
#        bloodfile=re.sub('.vcf','',syn_files[blood[0]])
#    else:
#        bloodfile=''

vcfdir='/mnt/huge/dermalWgs/varDict'
for outvcf in os.listdir(vcfdir):
        if 'vcf' not in outvcf or 'vep' in outvcf:
                continue
	outvcf=os.path.join(vcfdir,outvcf)
        cmfile=re.sub('.vcf','.sh',os.path.basename(outvcf))
        #patsh=open(cmfile,'w')
        ##first create vcf file with only two samples
#        bcftoolscmd="bcftools view %s -s %s -U"%(vcf_file,sstring)
#        if(bloodfile==''):
 #           outvcf="patient%s_tumor_%s_only.vcf"%(pat,t)
        outmaf=re.sub('.vcf','.maf',os.path.basename(outvcf))
        #annotationstr="'{\"dataType\":\"WGS\",\"tissueType\":\"tumorVsNormal\",\"patientId\":\""+tumfile_annotations['patientID'][0]+"\""
 #       annotationstr=annotationstr+",\"tissueID\":\""+tumfile_annotations['tissueID'][0]+"\"}'"


        ##now store paired VCF
       # synapse_upload_vcf="synapse store "+outvcf+".gz --parentId=syn5522791 --annotations "+annotationstr#+' --used '+usedstr

        #        patsh.write(synapse_upload_vcf+'\n\n')
        fdeets=re.sub(".vcf","",os.path.basename(outvcf)).split('_')
        pat='patient_'+fdeets[1]
        normid=pat+'_normal_'+fdeets[3]
        tumid='tumor_'+fdeets[5]
      #  patsh.write('bgzip -d '+outvcf+'.gz\n\n')
        vcf2maf_cmd=vcf2maf+" --input-vcf %s --vcf-tumor-id %s --tumor-id %s"%(outvcf,tumid,pat+'_'+tumid)
        vcf2maf_cmd+=' --vcf-normal-id %s --normal-id %s'%(normid,normid)
        vcf2maf_cmd+=" --output-maf %s --vep-forks 28 --species homo_sapiens --ref-fasta %s"%(outmaf,reffasta)
        print vcf2maf_cmd
	os.system(vcf2maf_cmd)
#        patsh.write(vcf2maf_cmd+'\ngzip '+outmaf+'\n')
        #then these file should be uploaded to synapse

#        synapse_upload_maf="synapse store "+outmaf+".gz --parentId=syn5522808 --annotations "+annotationstr#+' --used '+usedstr
#        patsh.write(synapse_upload_maf+'\n')

        #patsh.write(bcftoolscmd+'>'+outvcf+'\n'+vcf2maf_cmd+'\n')
        #os.system(cmd)
      #  patsh.close()
    #os.system('sh '+cmfile+' &')

        #examples for patient 1
'''
syn4987466 syn4985519
syn4985596 syn4985519
syn4985378 syn4985519
##here is the test command, as far as I can tell...
perl ../../../vcf2maf-master/vcf2maf.pl --input-vcf recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf --tumor-id=syn4987466 --normal-id=syn4985519 --output-maf=tumorVsNormal_pat1.maf
'''
