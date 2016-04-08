#!/usr/local/bin/python

'''
Goal is to download data from synapse and
format for input into comet2
'''
import os,sys,re
cometdir='~/comet/'

os.system('Rscript ./formatCometData.R')

afiles=[a for a in os.listdir('./') if 'tab' in a and 'mutGenes' in a]
heatfiles=[]
#for a in afiles:
#  altfile=re.sub('.tab','.heat',a)
#  print altfile
#  heatfiles.append(altfile)
#  os.system('python '+os.path.join(cometdir,'generateHeat.py')+' scores -o '+altfile+' -hf '+a)
heatfiles = afiles


for cfile in heatfiles:
    for set in ['2','3','4','5']:
        base=re.sub('.tab','SetsOf'+set,os.path.basename(cfile))
	
	cmdstr=' python ~/comet/run_comet_full.py -o %s --parallel -np 1000 -N 1000000 -m %s  -ks %s'%(base,cfile,set)
    	print cmdstr
    	os.system(cmdstr)
   
	#permute
	#now compute significance
	 
