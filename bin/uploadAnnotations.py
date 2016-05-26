import synapseclient
import synapseutils as synu 
import os
syn = synapseclient.login()

#Add in alternate ID to proteomics data
def addAlternateIdToProteomics():
	meta = syn.get('syn5361378')
	metadata = pd.read_csv(meta.path,sep="\t")

	DNAid = metadata['DNA Patient/tumor ID'] 
	RNAid = metadata['RNA Patient/tumor ID']
	RNAsample = metadata['RNA Sample ID']
	DNAsample = metadata['DNA Sample ID']

	DNAid = [i[-4:] for i in DNAid]
	RNAid = [i[-4:] for i in RNAid]

	fileList = syn.query('select id from file where parentId =="syn4984949" AND patientID != ""')

	for i in fileList['results']:
		temp = syn.get(i['file.id'], downloadFile=False)
		annots = temp.annotations
		index = DNAsample[DNAsample.str.contains(temp.annotations.sampleID[0])==True]
		alternate=""
		if len(index == 1):
			alternate = RNAid[index.index[0]]
		index = RNAsample[RNAsample.str.contains(temp.annotations.sampleID[0])==True]
		if len(index == 1):
			alternate = DNAid[index.index[0]]
		annots['alternateTumorID']=alternate
		syn.store(temp,forceVersion=False)


dermalNF = {"fundingAgency":"CTF",
 			"speciesName":"homoSapiens"}

def annotate_proteomics():
	###PROTEOMICS###
	directory = synu.walk(syn,"syn6101352")
	for dirpath, dirname, filename in directory:
		annots = {"topic":"proteomicsExperiment",
				  "disease":"NF1",
				  "diseaseSubtype":"dermalNF",
				  "tissueType":"neurofibroma",
				  "tissueSubtype":"dermalNeurofibroma",
				  "nf1Genotype":"-/-"}	
		if os.path.basename(dirpath[0]) == "Original Calls":
			annots.update({"dataType":"proteomicsData",
					  "dataSubtype":"peptideIdentification"
					  "fileFormat":"text",
					  "formatSubtype":"tsv"})
		else:
			annots.update({"dataType":"proteomicsData",
					  "dataSubtype":"proteinMatrix"
					  "fileFormat":"text",
					  "formatSubtype":"tsv"})
		proteomics = dict(dermalNF.items() + annots.items())
		for i in filename:
			temp = syn.get(i[1],downloadFile=False)
			temp.annotations.update(proteomics)
			for i in temp.annotations:
				if type(temp[i]) == list and len(temp[i]) == 1:
					temp[i] = temp[i][0] 
			syn.store(temp, forceVersion=False)

###RNA###
def annotate_rna():
	directory = synu.walk(syn,"syn6035832")
	for dirpath, dirname, filename in directory:
		annots = {"topic":"alignment",
				  "subTopic":"RNAseq"}
		if "ENCODE controls" in  dirpath[0]:
			annots.update({"disease":"normal",
					   "nf1Genotype":"+/+",
					   "tissueType":"skin"})
		else:
			annots.update({"disease":"NF1",
					   "diseaseSubtype":"dermalNF",
					   "tissueType":"neurofibroma",
					   "tissueSubtype":"dermalNeurofibroma",
					   "nf1Genotype":"-/-"})
		if os.path.basename(dirpath[0]) == "Cufflinks Quantitation":
			annots.update({"dataType":"geneExpressionData",
					   "dataSubtype":"geneCountsFile",
					   "fileType":"text"})
		elif os.path.basename(dirpath[0]) == "FeatureCounts Quantitation":
			annots.update({"dataType":"geneExpressionData",
					   "dataSubtype":"geneCountsFile",
					   "fileType":"text"})
		elif os.path.basename(dirpath[0]) == "Gene matrices":
			annots.update({"dataType":"geneExpressionData",
					   "dataSubtype":"geneExpressionMatrix",
					   "fileType":"text",
					   "fileSubtype":"csv"})
		elif os.path.basename(dirpath[0]) == "Dermal NF tumor alignments":
			annots.update({"dataType":"sequenceAlignment",
					   "dataSubtype":"rnaSequenceAlignment",
					   "fileType":"binary"})
		elif os.path.basename(dirpath[0]) == "BAM files":
			annots.update({"dataType":"sequenceAlignment",
					   "dataSubtype":"rnaSequenceAlignment",
					   "fileType":"binary"})
		elif os.path.basename(dirpath[0]) == "Pre-aligned quantitation":
			annots.update({"dataType":"geneExpressionData",
					   "dataSubtype":"geneCountsFile",
					   "fileType":"text",
					   "fileSubtype":"tsv"})
		final = dict(dermalNF.items() + annots.items())
		for i in filename:
			temp = syn.get(i[1],downloadFile=False)
			temp.annotations.update(final)
			if i[0].endswith(".gtf"):
				temp.fileSubtype = "gtf"
			elif i[0].endswith(".bam"):
				temp.fileSubtype = "BAM"
			elif i[0].endswith(".bai"):
				temp.fileSubtype = "BAI"
			elif temp.get('fileSubtype',None) is None:
				temp.fileSubtype = ''
			for i in temp.annotations:
				if type(temp[i]) == list and len(temp[i]) == 1:
					temp[i] = temp[i][0] 
			#print(temp.annotations)
			syn.store(temp, forceVersion=False)

###SNP###
def annotate_snp():
	directory = synu.walk(syn,"syn5004874")
	for dirpath, dirname, filename in directory:
		annots = {"topic":"geneticVariation",
				  "subTopic":"dnaPolymorphism",
				  "disease":"NF1",
				  "diseaseSubtype":"dermalNF",
				  "tissueType":"neurofibroma",
				  "tissueSubtype":"dermalNeurofibroma",
				  "nf1Genotype":"-/-"}
		final = dict(dermalNF.items() + annots.items())
		for i in filename:
			temp = syn.get(i[1],downloadFile=False)
			temp.annotations.update(final)
			if i[0].startswith("3096"):
				temp.dataType = "report"
				temp.dataSubtype = "snpArrayReport"
				temp.fileType = "text"
				temp.fileSubtype = "csv"
			else:
				temp.dataType = "annotation"
				temp.dataSubtype = "snpAnnotation"
				temp.fileType = "text"
				temp.fileSubtype = "csv"			
			for i in temp.annotations:
				if type(temp[i]) == list and len(temp[i]) == 1:
					temp[i] = temp[i][0] 
			#print(temp.annotations)
			syn.store(temp, forceVersion=False)

### Sample Info ###
def annotate_sampleInfo():
	directory = synu.walk(syn,"syn4984723")
	for dirpath, dirname, filename in directory:
		annots = {"topic":"dataIdentityAndMapping",
				  "subTopic":"clinicalVariables",
				  "disease":"NF1",
				  "diseaseSubtype":"dermalNF",
				  "tissueType":"neurofibroma",
				  "tissueSubtype":"dermalNeurofibroma",
				  "nf1Genotype":"-/-",
				  "dataType":"annotation",
				  "dataSubtype":"sampleAnnotation"}
		final = dict(dermalNF.items() + annots.items())
		for i in filename:
			temp = syn.get(i[1],downloadFile=False)
			temp.annotations.update(final)
			temp.fileSubtype = i[0].split(".")[-1]			
			for i in temp.annotations:
				if type(temp[i]) == list and len(temp[i]) == 1:
					temp[i] = temp[i][0] 
			#print(temp.annotations)
			syn.store(temp, forceVersion=False)

###WGS####
def annotate_WGS():
	directory = synu.walk(syn,"syn4984626")
	for dirpath, dirname, filename in directory:
		annots = {"topic":"geneticVariation",
				  "subTopic":"dnaMutation",
				  "disease":"NF1",
				  "diseaseSubtype":"dermalNF",
				  "tissueType":"neurofibroma",
				  "tissueSubtype":"dermalNeurofibroma",
				  "nf1Genotype":"-/-"}
		final = dict(dermalNF.items() + annots.items())
		for i in filename:
			temp = syn.get(i[1],downloadFile=False)
			temp.annotations.update(final)
			if i[0].endswith(".vcf.gz") or i[0].endswith(".vcf"):
				temp.fileFormat = "typedText"
				temp.formatSubtype = "VCF"
			elif i[0].endswith(".maf.gz") or i[0].endswith(".maf"):
				temp.fileFormat = "typedText"
				temp.formatSubtype = "MAF"				
			for i in temp.annotations:
				if type(temp[i]) == list and len(temp[i]) == 1:
					temp[i] = temp[i][0] 
			#print(temp.annotations)
			syn.store(temp, forceVersion=False)
