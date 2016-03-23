'''
Goal is to download files from synapse and iterate over different parameters
'''


import os,re
import synapseclient
import networkx

syn = synapseclient.Synapse()
syn.login()

allfiles=syn.query("SELECT name,id from entity where parentId =='syn5804586'")['results']

omics_int_dir='~/OmicsIntegrator/'

frac_ids=[a['entity.id'] for a in allfiles if 'fracOfProteins' in a['entity.name']]
mean_ids=[a['entity.id'] for a in allfiles if 'meanOfProteins' in a['entity.name']]
sum_ids=[a['entity.id'] for a in allfiles if 'sumOfProteins' in a['entity.name']]

#we're only evaluating mean values now.
allfiles={'proteinMean':mean_ids} #{'proteinFraction':frac_ids,'proteinMean':mean_ids,'proteinTotal':sum_ids}

def conf_prep(mu,beta,D,w):
    file = open("conf.txt","w")
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

mu_range = [0.009]
beta_range = [25]
w_range = [2]

edge_file = os.path.join(omics_int_dir,"data/iref_mitab_miscore_2013_08_12_interactome.txt")
conf_file = "conf.txt"
D = 5

outfolder='syn5816783'
this_script='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-22/runForestWithRand.py'

##lots of for loops to query everything involved, and then read in file to
##collect stats:
##1- individual number of edges#
##2- number of trees
##3- whther or not ubiquitin was included

##goal it to write it all into tab-delimited file so we can
##analyze it in R.

for a in allfiles.keys():
    output_path = "../2016-03-22/"+a
    if not os.path.exists(output_path):
        os.system('mkdir '+output_path)
    for pf in allfiles[a]:
        prize_file = syn.get(pf).path
        if 'AllSamples' in prize_file:
            prefix='proteinAcrossSamples'
        else:
            patient=re.sub('.tab','',os.path.basename(prize_file)).split('Patient')[1]
            prefix = 'proteinAcrossPatient%s'%(patient)
        for mu in mu_range:
            for beta in beta_range:
                for w in w_range:
                    out_label = "%s_w%f_beta%d_D%d_mu%f" %(prefix,w,beta,D,mu)
                    conf_prep(mu,beta,D,w)
                    out_label = "%s_w%f_beta%d_D%d_mu%f" %(prefix,w,beta,D,mu)
                    ne_cmd="python %s --prize %s --edge %s --conf conf.txt  --noisyEdges 100 --msgpath ~/msgsteiner-1.1/msgsteiner --outpath %s --outlabel %s" %(os.path.join(omics_int_dir,'scripts/forest.py'),prize_file,edge_file,output_path,out_label)
                    print ne_cmd
                    os.system(ne_cmd)

                    sp_cmd="python %s --prize %s --edge %s --conf conf.txt
                    --shuffledPrizes 100 --msgpath ~/msgsteiner-1.1/msgsteiner --outpath %s --outlabel %s"%(os.path.join(omics_int_dir,'scripts/forest.py'),prize_file,edge_file,output_path,out_label)
                    print sp_cmd
                    os.system(sp_cmd)
                    #opt_file=os.path.join(input_path,out_label+'_optimalForest.sif')
                    #dum_file=os.path.join(input_path,out_label+'_dummyForest.sif')
                    #if os.path.exists(opt_file) and os.path.exists(dum_file):
                    #    net_stats=netStats(opt_file,dum_file)
                    #    ofile.write("%s\t%s\t%f\t%d\t%f\t%d\t%d\t%d\t%s\t%s\n"%(a,prefix,mu,beta,w,net_stats['numEdges'],\
                     #                                                       net_stats['numTrees'],net_stats['numNodes'],\
                     #                                                       net_stats['treeSizes'],net_stats['hasUBC']))
                    outfiles=[]
                    for of in outfiles:
                        syn.store(synapseclient.File(of,parentId=outfoler),\
                                activityName='ran with randomization',used = [pf,{'url':this_script,'wasExecuted':True}])
