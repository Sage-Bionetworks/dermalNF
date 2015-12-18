import os,re,sys
import pandas
import math

first_bamdir='/scratch/DERMALNF/disk1'
second_bamdir='/scratch/DERMALNF/disk3/wgs_analysis/'

first_fis=[a for a in os.listdir(first_bamdir) if '.bam' in a]
second_fis=[a for a in os.listdir(second_bamdir) if '.bam' in a]

#now get tumor/normal pairs from table!
import synapseclient
syn = synapseclient.Synapse()
syn.login()

stcommand='samtools mpileup -f ../../../lib/ucsc.hg19.fasta '
vscommand='java -r ~/VarScan.v2.3.9.jar pileup2snp'

def getBamPath(synid):
    '''
    Quick function to find bamfile
    '''
    ns=df['WGS'][synid][0]
    f=syn.get(ns,downloadFile=False).name
    bamf=re.sub(".vcf",'.bam',f)
    if bamf in first_fis:
        return os.path.join(first_bamdir,bamf)
    elif bamf in second_fis:
        return os.path.join(second_bamdir,bamf)
    else:
        print '%s not found on scratch'%(bamf)
        return ''

allpats=['1','2','3','4','5','6','7','8','9','10','11','12','13']
for p in allpats:
    res=syn.tableQuery("SELECT DnaID,WGS FROM syn5556216 where Patient=%s"%(p))
    normind=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if a=='PBMC']
    if len(normind)==0:
        print 'No normal file found for patient '+p
        next
    normfile=getBamPath(normind)

    tuminds=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if (a!='PBMC' and not math.isnan(float(a)))]
    for tu in tums:
        #run pileup command
        bf=getBamPath(tu)
        outfile='patient_%s_tumor_%s_vs_normal_%s.snp'%(p,tu,normind)
        print outfile
        fullcmd=stcommand+bf+' '+normfile+' | '+vscommand+'>'+outfile
        print fullcmd
