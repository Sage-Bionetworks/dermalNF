'''
Try to rename headeres on merged VCF file
'''

import synapseclient,re,os

syn=synapseclient.Synapse()
syn.login()

#query for original patient/sample identifiers
res=syn.query("select id,name from entity where parentId=='syn4984931'")

vcffiles=dict([(r['entity.name'],r['entity.id']) for r in res['results'] if 'hard-filtered' not in r['entity.name']])

sampleids={}
for f in vcffiles.keys():
    si=vcffiles[f]
    annotes=syn.getAnnotations(si)
    newname='patient_'+re.sub('CT0+','',annotes['patientID'][0])+'_'+annotes['tissueType'][0]+'_'+annotes['tissueID'][0]
    print newname
    sampleids[re.sub('.vcf','',f)]=newname

##now create mapping from file name
##need to get actual column name and header from VCF
olheader=open('../2015-11-24/headerFile.txt','r').readlines()
olheader[len(olheader)-1]

sampstring=olheader[len(olheader)-1].strip().split('\t')
newstring=sampstring

for i in range(9,len(sampstring)):
    slnum=sampstring[i]
    ##match sample number to patient and tumor
    if slnum in sampleids.keys():
        newnum=sampleids[slnum]
    else:
        newnum=slnum
    newstring[i]=newnum

newheader=olheader
newheader[len(olheader)-1]='\t'.join(newstring)+'\n'

open('newHeaderFile.txt','w').writelines(newheader)

#now we get to reheader!!!!
newcmd='bcftools reheader ../2015-11-24/geneRegionsWH.vcf -h newHeaderFile.txt -o geneRegionsOnlyWithSampleNames.vcf'
print newcmd
os.system(newcmd)
