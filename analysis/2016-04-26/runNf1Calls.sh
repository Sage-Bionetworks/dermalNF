#python ../../bin/varDictVcfDataProcessing.py /mnt/huge/dermalWgs/vdj

cmd='perl /home/ubuntu/VarDictJava/VarDict/vcf2txt.pl'
files='/mnt/huge/dermalWgs/vdj/*[0-9].snpeff.vcf'
$cmd $files > nf1_combined.maf

#for file in /mnt/huge/dermalWgs/vdj/*vcf
#do
#java -jar ~/

#$cmd $file >> nf1_all_merged.maf
#done

