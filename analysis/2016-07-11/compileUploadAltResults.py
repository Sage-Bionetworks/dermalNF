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

nf1reg='-R chr17:29421944-29704695:NF1'

vd_dir='/home/ubuntu/VarDictJava/'
def runVarDictOnBed(normfile,tumfile,normsamp,tumsamp,bedfile,suffix='',addPre=True,cmdfile=''):
    reference='/home/ubuntu/dermalNF/lib/ucsc.hg19.fasta'
    #<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | VarDict/var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
   # bedfile = nf1reg
    vdex=os.path.join(vd_dir,'build/install/VarDict/bin/VarDict')
    vdcmd="vdcmd=\"%s -G %s -f 0.01 -N %s -b \\\"%s|%s\\\" -c 1 -S 2 -E 3"%(vdex,reference,tumsamp,tumfile,normfile)
    tscmd="tscmd=\"%s/testsomatic.R\""%(os.path.join(vd_dir,'VarDict'))
    vcfcmd="vcfcmd=\"%s/var2vcf_paired.pl -N \\\"%s|%s\\\" -f 0.01\""%(os.path.join(vd_dir,'VarDict'),tumsamp,normsamp)
    outpre=normsamp+'_'+tumsamp+'.vcf'+suffix
    if cmdfile=='':
        if addPre:
	    os.system(tscmd+';'+vcfcmd)
        newvd=vdcmd+' -g 4 %s\"'%(bedfile)
        os.system(newvd)
        os.system('$vdcmd|$tscmd|$vcfcmd>'+outpre)
    else:
	if addPre:
             cmdfile.write(tscmd+'\n'+vcfcmd+'\n')
        newvd=vdcmd+' -g 4 %s\"'%(bedfile)
        cmdfile.write(newvd+'\n')
        cmdfile.write('$vdcmd|$tscmd|$vcfcmd>%s\n'%(outpre))
    return outpre

def runSnpEff(vcf,cmdfile=''):
    cmd='java -jar /home/ubuntu/snpEff/snpEff.jar hg19 %s'%(vcf)
    newvcf=re.sub('vcf','snpeff.vcf',vcf)
    if cmdfile=='':
	os.system(cmd+'>'+newvcf)
    else:
	cmdfile.write(cmd+'>'+newvcf+'\n')

    return newvcf


def runVarDict(normfile,tumfile,normsamp,tumsamp,cmdfile=''):
    reference='/home/ubuntu/dermalNF/lib/ucsc.hg19.fasta'

    vdcmd="vdcmd=\"/home/ubuntu/VarDict/vardict -G %s -f 0.01 -N %s -b %s|%s -c 1 -S 2 -E 3"%(reference,normsamp,normfile,tumfile)
    tscmd="tscmd=\"/home/ubuntu/VarDict/testsomatic.R\""
    vcfcmd="vcfcmd=\"/home/ubuntu/VarDict/var2vcf_paired.pl -N %s|%s -f 0.01\""%(normsamp,tumsamp)
    outpre=normsamp+'_'+tumsamp+'.vcf'
    bf=range(1,21)
    if cmdfile=='':
        os.system(tscmd+';'+vcfcmd)
    else:
	cmdfile.write(tscmd+'\n'+vcfcmd+'\n')

    for b in bf:
	bedfile='/home/ubuntu/VarDict/hg19_knownGeneRed.bed.%d'%(b)
        vcf=runVarDictOnBed(normfile,tumfile,normsamp,tumsamp,bedfile,suffix='.%d'%(b),addPre=False,cmdfile=cmdfile)
        newvcf=runSnpEff(vcf,cmdfile)
#	allvcf=' '.join([outpre+'.'+str(b) for b in bf])
#        os.system('bcftools merge '+allvcf+' -o '+outpre)

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


def compileVcfs(base,numFiles,cmdfile=''):
    outfiles=[base+'.%d'%(i) for i in range(1,numFiles+1)]
    base=os.path.basename(base)
    for o in outfiles:
	if not os.path.exists(o+'.gz'):
 	    os.system('bgzip '+o)
    	os.system('bcftools index -f '+o+'.gz')
    cmd='bcftools concat -a '+' '.join([o+'.gz' for o in outfiles])+' -o '+base
    print cmd
    if cmdfile=='':
	if not os.path.exists(base):
            os.system(cmd)
	if not os.path.exists(base+'.gz'):
     	    os.system('bgzip '+base)
    else:
        cmdfile.write(cmdfile+'\n')
    return base+'.gz'

def applyVcfFilter(vcf,cmdfile=''):
    ofile=re.sub('vcf','filtered.vcf',vcf)
    cmd='bgzip '+vcf+'\nbcftools index -t '+vcf+'.gz\nbcftools view '+vcf+'.gz -o '+ofile+'.gz -f.,PASS -Oz'
    print cmd
    os.system(cmd)
    return ofile+'.gz'

def makeMafFromVcf(vcffile,cmdfile=''):
    cmd='perl ~/VarDictJava/VarDict/vcf2txt.pl'
    maffile=re.sub(".vcf.gz",".maf",vcffile)
    cmd=cmd+' '+vcffile
    print cmd
    if cmdfile=='':
        os.system(cmd+'>'+maffile)
	os.system('bgzip '+maffile)
    else:
        cmdfile.write(cmdfile+'>'+maffile+'\n')

    return maffile+'.gz'

allpats=['5','6','8','9']
#allpats = ['1']
for p in allpats:
    res=syn.tableQuery("SELECT DnaID,WGS FROM syn5556216 where Patient=%s"%(p))
    df=res.asDataFrame()
    normind=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if a=='PBMC']
    # if len(normind)==0:
    #     print 'No normal file found for patient '+p
    #     next
    # normscript=open('patient'+p+'_'+normind[0]+'bamProcessing.sh','w')
    # normfile=updateBams(getBamPath(normind[0],normscript),normscript)
    # normscript.write('cat "\\npatient '+p+' complete">> normstats\n')
    # normscript.close()

    tuminds=[df['WGS'][ind] for ind,a in enumerate(df['DnaID']) if (a!='PBMC' and not math.isnan(float(a)))]
#    tuminds=tuminds[0:1]
#    print tuminds
    for tu in tuminds:
        if isinstance(tu,float):
            continue
      #run pileup command
 #       cmdfile=open(tu+'_patient_upload_'+p+'_cmd.sh','w')
 #       cmdfile.write('#/bin/bash\n')
        cmdfile=''
#        bf=updateBams(getBamPath(tu,cmdfile),cmdfile)
#        outfile='patient_%s_tumor_%s_vs_normal_%s.snp'%(p,tu,normind[0])
	#if tu not in ["syn4985400","syn4985585","syn4984995","syn4985566","syn4985454","syn4985458","syn4984997"]:
	 #   print tu+' not in selected list, cmoving on'
	 #   continue
        normsamp='normal_%s'%(normind[0])
        tumsamp='patient_%s_tumor_%s'%(p,tu)

        base='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/'
        script=base+'2016-07-11/compileUploadAltResults.py'

       	svcf ='/mnt/huge/dermalWgs/varScan/'+tumsamp+'_vs_'+normsamp+'_VarScan_somatic.snp.snpeff.vcf'
        ivcf = '/mnt/huge/dermalWgs/varScan/'+tumsamp+'_vs_'+normsamp+'_VarScan_somatic.indel.snpeff.vcf'


    for vcf in [svcf, ivcf]:
        #first do snp file
        sv_file = synapseclient.File(vcf,description='Tumor-normal SNP VCF file',parentId='syn6834340')
        sv_act = synapseclient.activity.Activity(name='Vardict and snpeff',
                                                 description = 'Recalculated variants from bam files for Chromosome 17',
                                                used=[normind[0],tu],executed=script)
        si=syn.store(sv_file,activity=sv_act)


        fvcf = applyVcfFilter(vcf)

        sf_file = synapseclient.File(fvcf,description='Tumor-normal Filtered VCF file',parentId='syn6834340')
        sf_act=synapseclient.activity.Activity(name='Filtered vcf from varscan and snpeff',used=[normind[0],tu,si],executed=script)
        si=syn.store(sf_file,activity=sf_act)

        maffile = makeMafFromVcf(fvcf,cmdfile=cmdfile)

        sm_file = synapseclient.File(maffile,description='Tumor-normal filtered MAF file',parentId='syn6834373')
        sm_act=synapseclient.activity.Activity(name='MAF from varscan and snpeff',used=[normind[0],tu,si],executed=script)
        syn.store(sm_file,activity=sm_act)
