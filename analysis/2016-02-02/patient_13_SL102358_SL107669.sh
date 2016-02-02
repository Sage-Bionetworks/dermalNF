bcftools view /home/ubuntu/dermalNF/analysis/2015-12-16/recalibrated_indels_for_recalibrated_snps_for_dermalNFmerged.vcf.gz -s SL102358,SL107669 -U>patient13_tumor_syn4985372_vs_norm_syn4985529.vcf
bgzip patient13_tumor_syn4985372_vs_norm_syn4985529.vcf
synapse store patient13_tumor_syn4985372_vs_norm_syn4985529.vcf.gz --parentId=syn5522791 --annotations '{"dataType":"WGS","tissueType":"tumorVsNormal","patientId":"CT000013","tissueID":"0007"}'

bgzip -d patient13_tumor_syn4985372_vs_norm_syn4985529.vcf.gz

perl ../../../vcf2maf-master/vcf2maf.pl --input-vcf patient13_tumor_syn4985372_vs_norm_syn4985529.vcf --vcf-tumor-id SL102358 --vcf-tumor-id SL107669 --output-maf tumorVsNormal_pat13_tumor_syn4985372_vs_norm_syn4985529.maf --vep-forks 9 --species homo_sapiens --ref-fasta ../../lib/ucsc.hg19.fasta
gzip tumorVsNormal_pat13_tumor_syn4985372_vs_norm_syn4985529.maf
synapse store tumorVsNormal_pat13_tumor_syn4985372_vs_norm_syn4985529.maf.gz --parentId=syn5522808 --annotations '{"dataType":"WGS","tissueType":"tumorVsNormal","patientId":"CT000013","tissueID":"0007"}'

