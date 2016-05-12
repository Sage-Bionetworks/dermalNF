import os,re,sys
import pandas
import math


bamrepo='s3://ctf-data/wgs-bams/'

#now get tumor/normal pairs from table!
import synapseclient
syn = synapseclient.Synapse()
syn.login()

stcommand='samtools mpileup -f ../../../lib/ucsc.hg19.fasta '
vscommand='java -r ~/VarScan.v2.3.9.jar pileup2snp'

def runMutect(normfile,tumfile,out_prefix):
    reference='../../lib/ucsc.hg19.fasta'
    cosmic='../../lib/sorted_b37_cosmic_v54_120711_withchr.vcf'
    dbsnp='../../lib/sorted_dbsnp_138.hg19.vcf'
    intervals='../../lib/ucsc.hg19.dict'
    cmd='java -jar ~/muTect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence %s --cosmic %s \
         --dbsnp %s --input_file:normal %s --input_file:tumor %s --out %s -cov %s'%(reference,cosmic,dbsnp,normfile,tumfile,out_prefix+'.out',out_prefix+'_coverage.wig.txt')
    print cmd


def updateBams(bamfile): 
    outfile=re.sub('.bam','_rg.bam',bamfile)
    picmd='java -jar ~/picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups \
           I=%s O=%s RGID=1 RGLB=hg19 RGPL=illumina RGPU=dragen RGSM=%s'%(bamfile,outfile,re.sub('.bam','',bamfile))
    print picmd
    if not os.path.exists(outfile):
	os.system(picmd)

    ordered=re.sub('.bam','_ordered.bam',outfile)
    pic2cmd='java -jar ~/picard-tools-2.1.1/picard.jar ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s'%(outfile,ordered,'../../lib/ucsc.hg19.fasta')
    if not os.path.exists(ordered):
	print pic2cmd
	os.system(pic2cmd)

    return ordered


def getBamPath(synid):
    '''
    Quick function to find bamfile
    '''
    f=syn.get(synid,downloadFile=False).name
    bamf=re.sub(".vcf",'.bam',f)
    awscmd='aws s3 cp %s %s'%(os.path.join(bamrepo,bamf),bamf)
    print awscmd
    if not os.path.exists(bamf):
    	os.system(awscmd)
    
    return bamf

allpats=['1','2','3','4','5','6','8','9','11']
allpats = ['1']
for p in allpats:
    res=syn.tableQuery("SELECT DnaID,WGS FROM syn5556216 where Patient=%s"%(p))
    df=res.asDataFrame()
    normind=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if a=='PBMC']
    if len(normind)==0:
        print 'No normal file found for patient '+p
        next
    normfile=updateBams(getBamPath(normind[0]))

    tuminds=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if (a!='PBMC' and not math.isnan(float(a)))]
    tuminds=tuminds[0:1]
    for tu in tuminds:
        #run pileup command
        bf=updateBams(getBamPath(tu))
        outfile='patient_%s_tumor_%s_vs_normal_%s.snp'%(p,tu,normind[0])
        print outfile
        cmd = runMutect(normfile,bf,outfile)
