##first we need to move VCF files into their own directory, and upload
##the merged VCF file

import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()

wgs_vcf='syn4984931'

query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")
##now get metadata for all files
syn_files=dict([(res['entity.id'],res['entity.name']) for res in \
                query_res['results'] if 'hard-filtered' not in \
                res['entity.name'] and '.vcf' in res['entity.name']])

syn_ids=syn_files.keys()
#print syn_files

uf=syn_ids+['https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/wgsVariantRefinement.py']
#lastly, add in GATK processing file
gafile=synapseclient.entity.File('../2015-11-20/recalibrated_indels_for_recalibrated_snps_for_dermalNFall.vcf.gz',parentId='syn5522790',used=uf)
syn.store(gafile)

npid='syn5522788'
for sy in syn_ids:
    sf=syn.get(sy,downloadFile=False)
    sf.properties.parentId=npid
    syn.store(sf)



#now get all hard-filtered
syn_files=dict([(res['entity.id'],res['entity.name']) for res in \
                query_res['results'] if 'hard-filtered' in \
                res['entity.name'] and '.vcf' in res['entity.name']])

syn_ids=syn_files.keys()
for sy in syn_ids:
    syn.delete(sy)
