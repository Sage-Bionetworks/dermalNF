import os,re,sys
import pandas
import math


bamrepo='s3://ctf-data/wgs-bams/'

#now get tumor/normal pairs from table!
import synapseclient
syn = synapseclient.Synapse()
syn.login()

#stcommand='samtools mpileup -f ../../../lib/ucsc.hg19.fasta '
#vscommand='java -r ~/VarScan.v2.3.9.jar pileup2snp'

def runVarDict(normfile,tumfile,normsamp,tumsamp,cmdfile=''):
    reference='/home/ubuntu/dermalNF/lib/ucsc.hg19.fasta'

    vdcmd="vdcmd=\"/home/ubuntu/VarDict/vardict -G %s -f 0.01 -N %s -b %s|%s -c 1 -S 2 -E 3"%(reference,normsamp,normfile,tumfile)
    tscmd="tscmd=\"/home/ubuntu/VarDict/testsomatic.R\""
    vcfcmd="vcfcmd=\"/home/ubuntu/VarDict/var2vcf_somatic.pl -N %s|%s -f 0.01\""%(normsamp,tumsamp)
    outpre=normsamp+'_'+tumsamp+'.vcf'
    bf=range(1,11)
    if cmdfile=='':
        os.system(tscmd+';'+vcfcmd)
        for b in bf:
	    bedfile='/home/ubuntu/VarDict/hg19_knownGene.bed.%d'%(b)    
            newvd=vdcmd+' -g 3 %s\"'%(bedfile)
            os.system(newvd)
            os.system('$vdcmd|$tscmd|$vcfcmd>'+outpre+'.'+str(b))
	allvcf=' '.join([outpre+'.'+str(b) for b in bf])
        os.system('bcftools merge '+allvcf+' -o '+outpre)
    else:
        cmdfile.write(tscmd+'\n'+vcfcmd+'\n')
        for b in bf:
	    bedfile='/home/ubuntu/VarDict/hg19_knownGene.bed.%d'%(b)    
            newvd=vdcmd+' -g 3 %s\"'%(bedfile)
            cmdfile.write(newvd+'\n')
            cmdfile.write('$vdcmd|$tscmd|$vcfcmd>%s.%d\n'%(outpre,b))
	allvcf=' '.join([outpre+'.'+str(b) for b in bf])
        cmdfile.write('bcftools merge '+allvcf+' -o '+outpre+'\n')

def runMutect(normfile,tumfile,out_prefix,cmdfile=''):
    reference='/home/ubuntu/dermalNF/lib/ucsc.hg19.fasta'
    cosmic='/home/ubuntu/dermalNF/lib/sorted_b37_cosmic_v54_120711_withchr.vcf'
    dbsnp='/home/ubuntu/dermalNF/lib/sorted_dbsnp_138.hg19.vcf'
    intervals=['chr1','chr2','chr3','chr4','chr6','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']#'../../lib/ucsc.hg19.dict'
    if cmdfile=='':
	os.system('sudo apt-get -f install oracle-java7-set-default')
    else:
	cmdfile.write('sudo apt-get -f install oracle-java7-set-default\n')
    for iv in intervals:
	cmd='java -jar ~/GenomeAnalysisTK.jar -T MuTect2 -R %s -I:normal %s -I:tumor %s --out %s -L %s --cosmic %s --dbsnp %s'%(reference,normfile,tumfile,out_prefix+'_'+iv+'.vcf',iv,cosmic,dbsnp)

   # cmd='java -jar ~/muTect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence %s --cosmic %s \
   #      --dbsnp %s --input_file:normal %s --input_file:tumor %s --out %s -cov %s'%(reference,cosmic,dbsnp,normfile,tumfile,out_prefix+'.out',out_prefix+'_coverage.wig.txt')
    	pstr='Running muTectv2 on %s and %s on region %s'%(normfile,tumfile,iv)
   	if cmdfile=='':
       	    print pstr
      	    os.system(cmd)

    	else:
            cmdfile.write('echo "'+pstr+'"\n'+cmd+'\n')
    allfiles=' '.join([out_prefix+'_'+chr+'.vcf' for chr in intervals])
    allidx=' '.join([out_prefix+'_'+chr+'.vcf.idx' for chr in intervals])
    if cmdfile=='':
	os.system('bcftools merge '+allfiles+' -o '+out_prefix+'.vcf')
	os.system('rm '+allfiles)
	os.system('rm '+allidx)
    else:
	cmdfile.write('bcftools merge '+allfiles+' -o '+out_prefix+'.vcf\nrm '+allfiles+'\nrm '+allidx)


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
	    cmdfile.write('echo "Adding read groups to make %s"\n%s\n%s\n'%(outfile,picmd,'rm '+bamfile))

    pic2cmd='java -jar /home/ubuntu/picard-tools-2.1.1/picard.jar ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s TMP_DIR=/mnt/huge/tmp'%(outfile,ordered,'/home/ubuntu/dermalNF/lib/ucsc.hg19.fasta')
    if not os.path.exists(ordered):
	ustr='Re-ordering BAM file to make %s'%(ordered)
	if cmdfile=='':
	    print ustr
	    os.system(pic2cmd)
	    #ios.system('rm '+outfile)
	else:
	    cmdfile.write('echo "'+ustr+'"\n'+pic2cmd+'\n')#rm '+outfile+'\n')


    ind='samtools index %s'%(ordered)
    ustr='Running samtools to index %s'%(ordered)
    if not os.path.exists(ordered+'.bai'):
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
    bamf=re.sub(".vcf",'.bam',f)
    bampath='../'+bamf
    awscmd='aws s3 cp %s %s'%(os.path.join(bamrepo,bamf),bampath)
    #awscmd='aws s3 sync %s . --include %s'%(bamrepo,bamf)
 #   print awscmd
    if not os.path.exists(bampath) and not os.path.exists(re.sub('.bam','_rg_ordered.bam',bampath)) and not os.path.exists(re.sub('.bam','_rg.bam',bampath)):
	if cmdfile=='':
    	    os.system(awscmd)
	else:
	    cmdfile.write('echo "Getting %s from AWS"\n%s\n'%(bamf,awscmd))

    return bampath

allpats=['1','2','3','4','5','6','8','9','11']
for p in allpats:
    res=syn.tableQuery("SELECT DnaID,WGS FROM syn5556216 where Patient=%s"%(p))
    df=res.asDataFrame()
    normind=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if a=='PBMC']
    if len(normind)==0:
        print 'No normal file found for patient '+p
        next
    normscript=open('patient'+p+'_'+normind[0]+'bamProcessing.sh','w')
    normfile=updateBams(getBamPath(normind[0],normscript),normscript)
    normscript.write('cat "\\npatient '+p+' complete">> normstats\n')
    normscript.close()



    tuminds=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if (a!='PBMC' and not math.isnan(float(a)))]
#    tuminds=tuminds[0:1]
    print tuminds
    for tu in tuminds:
	if isinstance(tu,float):
	    continue
      #run pileup command
	cmdfile=open(tu+'_patient'+p+'_cmd.sh','w')
        cmdfile.write('#/bin/bash\n')
        bf=updateBams(getBamPath(tu,cmdfile),cmdfile)
        outfile='patient_%s_tumor_%s_vs_normal_%s.snp'%(p,tu,normind[0])
        normsamp='patient_%s_normal_%s'%(p,normind[0])
        tumsamp='tumor_%s'%(tu)

        cmd = runVarDict(normfile,bf,normsamp,tumsamp,cmdfile)
	cmdfile.close()
