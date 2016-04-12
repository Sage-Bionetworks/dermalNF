#!/usr/local/bin/python

'''
Goal is to download data from synapse and
format for input into hotnet2
'''
import os,sys,re
hotnetdir='~/hotnet2/'

os.system('Rscript ../../bin/formatHotnetData.R')

afiles=[a for a in os.listdir('./') if 'tab' in a]
heatfiles=[]
#for a in afiles:
#  altfile=re.sub('.tab','.heat',a)
#  print altfile
#  heatfiles.append(altfile)
#  os.system('python '+os.path.join(hotnetdir,'generateHeat.py')+' scores -o '+altfile+' -hf '+a)
heatfiles = afiles


for cfile in heatfiles:
    base='hprd_'+re.sub('.tab','',os.path.basename(cfile))
    cmdstr="python %s -hf %s -if %s -mf %s -pnp %s -ef %s -nn %s -o %s -c 2"%(os.path.join(hotnetdir,'runHotNet2.py'),cfile,os.path.join(hotnetdir,'influence_matrices/hprd/hprd_index_genes'),os.path.join(hotnetdir,"influence_matrices/hprd/hprd_ppr_0.4.h5"),os.path.join(hotnetdir, "influence_matrices/hprd/permuted/##NUM##/hprd_ppr_0.4.h5"),os.path.join(hotnetdir,'influence_matrices/hprd/hprd_edge_list'),base,base)
    print cmdstr
    os.system(cmdstr)
    
    base='iref_'+re.sub('.tab','',os.path.basename(cfile))
    cmdstr="python %s -hf %s -if %s -mf %s -pnp %s -ef %s -nn %s -o %s -c 2"%(os.path.join(hotnetdir,'runHotNet2.py'),cfile,os.path.join(hotnetdir,'influence_matrices/irefindex/iref_index_genes'),os.path.join(hotnetdir,"influence_matrices/irefindex/iref_ppr_0.45.h5"),os.path.join(hotnetdir, "influence_matrices/irefindex/permuted/##NUM##/iref_ppr_0.45.h5"),os.path.join(hotnetdir,'influence_matrices/irefindex/iref_edge_list'),base,base)
    print cmdstr
    os.system(cmdstr)
