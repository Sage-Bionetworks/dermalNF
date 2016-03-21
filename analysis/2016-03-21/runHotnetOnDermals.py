#!/usr/local/bin/python

'''
Goal is to download data from synapse and
format for input into hotnet2
'''
import os,sys,re
hotnetdir='../../../../src/hotnet2-1.0.1/'

os.system('Rscript ./formatHotnetData.R')

afiles=[a for a in os.listdir('./') if 'tab' in a]
heatfiles=[]
for a in afiles:
  altfile=re.sub('.tab','.heat',a)
  print altfile
  heatfiles.append(altfile)
  os.system('python '+os.path.join(hotnetdir,'generateHeat.py')+' scores -o '+altfile+' -hf '+a)
 # os.system('python '+os.path.join(hotnetdir,'generateHeat.py')+' mutation '+a)


for cfile in heatfiles:
    cmdstr="python %s -hf %s -if %s -mf %s -pnp %s -ef %s  "%(os.path.join(hotnetdir,'runHotNet2.py'),cfile,os.path.join(hotnetdir,'influence_matrices/irefindex/iref_index_genes'),os.path.join(hotnetdir,"influence_matrices/irefindex/iref_ppr_0.55.h5"),os.path.join(hotnetdir, "influence_matrices/irefindex/iref_edgelist_##NUM##"),os.path.join(hotnetdir,'influence_matrices/irefindex/iref_edge_list'))
    print cmdstr
    os.system(cmdstr)
