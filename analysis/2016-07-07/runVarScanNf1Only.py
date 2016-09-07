import os,re,sys
import pandas
import math


bamrepo='s3://ctf-data/wgs-bams/'

#now get tumor/normal pairs from table!
import synapseclient
syn = synapseclient.Synapse()
syn.login()

stcommand='samtools mpileup -f ~/dermalNF/lib/ucsc.hg19.fasta -l ./UCSC_hg19_knwonGene_2015_11_14_chr17.bed '
vscommand='java -r ~/VarScan.v2.3.9.jar somatic'


def runSnpEff(vcf,cmdfile=''):
    cmd='java -jar /home/ubuntu/snpEff/snpEff.jar hg19 %s'%(vcf)
    newvcf=re.sub('vcf','snpeff.vcf',vcf)
    if cmdfile=='':
        os.system(cmd+'>'+newvcf)
    else:
        cmdfile.write(cmd+'>'+newvcf+'\n')
    return newvcf


def runMutect(normfile,tumfile,out_prefix,cmdfile=''):
    reference='~/dermalNF/lib/ucsc.hg19.fasta'
    cosmic='~/dermalNF/lib/sorted_b37_cosmic_v54_120711_withchr.vcf'
    dbsnp='~/dermalNF/lib/sorted_dbsnp_138.hg19.vcf'

#intervals=['chr1','chr2','chr3','chr4','chr6','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']#'../../lib/ucsc.hg19.dict'
    if cmdfile=='':
        os.system('sudo apt-get -f install oracle-java7-set-default')
    else:
        cmdfile.write('sudo apt-get -f install oracle-java7-set-default\n')
 #   for iv in intervals:
    iv='chr17'
    vcffile=out_prefix+'_'+iv+'.vcf'
    cmd='java -jar ~/GenomeAnalysisTK.jar --analysis_type MuTect2 --reference_sequence %s --input_file:normal %s --input_file:tumor %s --out %s -L %s'%(reference,normfile,tumfile,vcffile,iv)

    pstr='Running muTectv2 on %s and %s on region %s'%(normfile,tumfile,iv)
    if not os.path.exists(vcffile):
        if cmdfile=='':
            print pstr
            os.system(cmd)
        else:
            cmdfile.write('echo "'+pstr+'"\n'+cmd+'\n')
    return vcffile


def updateBams(bamfile,cmdfile=''):
    '''
    The BAM files require adding read groups, re-ordering, and indexing before they can be run by MuTect
    '''
    outfile=re.sub('.bam','_rg.bam',bamfile)
    ordered=re.sub('.bam','_ordered.bam',outfile)
    picmd='sudo apt-get -f install oracle-java8-set-default\njava -jar ~/picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I=%s O=%s RGID=1 RGLB=hg19 RGPL=illumina RGPU=dragen RGSM=%s'%(bamfile,outfile,re.sub('.bam','',bamfile))
#    print picmd
    if not os.path.exists(outfile) and not os.path.exists(ordered):
    	if cmdfile=='':
	    print 'Adding read groups to make %s'%(outfile)
	    os.system(picmd)
	    os.system('rm '+bamfile)
	else:
	    cmdfile.write('echo "Adding read groups to make %s"\n%s\n'%(outfile,picmd))

    pic2cmd='java -jar ~/picard-tools-2.1.1/picard.jar ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s'%(outfile,ordered,'../../lib/ucsc.hg19.fasta')
    if not os.path.exists(ordered):
	ustr='Re-ordering BAM file to make %s'%(ordered)
	if cmdfile=='':
	    print ustr
	    os.system(pic2cmd)
	    os.system('rm '+outfile)
	else:
	    cmdfile.write('echo "'+ustr+'"\n'+pic2cmd+'\nrm '+outfile+'\n')


    ind='samtools index %s'%(ordered)
    if not os.path.exists(ordered+'.bai'):
    	ustr='Running samtools to index %s'%(ordered)
   	if cmdfile=='':
            print ustr
            os.system(ind)
    	else:
            cmdfile.write('echo "'+ustr+'"\n'+ind+'\n')

    return ordered


def getBamPath(synid,cmdfile=''):
    '''
    Quick function to find bamfile
    '''
    f=syn.get(synid,downloadFile=False).name
    bamf=os.path.join('../',re.sub(".vcf",'.bam',f))
    awscmd='aws s3 cp %s %s'%(os.path.join(bamrepo,bamf),bamf)

    outfile=re.sub('.bam','_rg.bam',bamf)
    ordered=re.sub('.bam','_ordered.bam',outfile)
    if not os.path.exists(bamf) and not os.path.exists(outfile) and not os.path.exists(ordered):
	if cmdfile=='':
	    os.system(awscmd)
	else:
	    cmdfile.write('echo "Getting %s from AWS"\n%s\n'%(bamf,awscmd))

    return bamf

allpats=['1','2','3','4','5','6','8','9','11']
for p in allpats:
    res=syn.tableQuery("SELECT DnaID,WGS FROM syn5556216 where Patient=%s"%(p))
    df=res.asDataFrame()
    normind=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if a=='PBMC']
    if len(normind)==0:
        print 'No normal file found for patient '+p
        next
    cmdfile=open('patient_'+p+'_cmds.sh','w')
    cmdfile.write('#/bin/bash\n')

    normfile=updateBams(getBamPath(normind[0],cmdfile),cmdfile)
    normpu='patient_%s_normal.pileup'%(p)
    stnorm=stcommand+normfile+'>'+normpu
    if not os.path.exists(normpu):
        cmdfile.write(stnorm+'\n')

    tuminds=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if (a!='PBMC' and not math.isnan(float(a)))]
    for tu in tuminds:
        #run pileup command
        if isinstance(tu,float):
            continue
        #run pileup command
        bf=updateBams(getBamPath(tu,cmdfile),cmdfile)
        outfile='patient_%s_tumor_%s_vs_normal_%s_VarScan_somatic'%(p,tu,normind[0])
        tumpu='patient_%s_tumor_%s.pileup'%(p,tu)
        #    print outfile
        #cmd = runMutect(normfile,bf,outfile,cmdfile)
        #newvcf = runSnpEff(cmd,cmdfile)

#        bf=getBamPath(tu)
#       outfile='patient_%s_tumor_%s_vs_normal_%s.snp'%(p,tu,normind)
        print outfile

        sttum=stcommand+bf+'>'+tumpu
        vstum=vscommand+' '+normpu+' '+tumpu+' '+outfile
#        print fullcmd
        cmdfile.write(sttum+'\n'+vstum+'\n')
    cmdfile.close()
