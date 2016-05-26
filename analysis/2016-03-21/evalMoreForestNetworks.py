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

allfiles={'proteinFraction':frac_ids,'proteinMean':mean_ids,'proteinTotal':sum_ids}

def conf_prep(mu,beta,D,w):
    file = open("conf.txt","w")
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

mu_range = [0,0.009,0.015]
beta_range = [1,5,10,15,20,25]
w_range = [1,2,3]

edge_file = os.path.join(omics_int_dir,"data/iref_mitab_miscore_2013_08_12_interactome.txt")
conf_file = "conf.txt"
D = 5

'''
Get network statistics for a particular .sif file to try to determine best
one to use!
'''
def netStats(opt_fname,dummy_fname):
    print 'Evaluating statistics from %s,%s'%(opt_fname,dummy_fname)

    opt_tab=[a.strip().split('\tpp\t') for a in open(opt_fname,'r').readlines()]
    dummy_tab=[[a.strip().split('\t')[0],a.strip().split('\t')[2]] for a in open(dummy_fname,'r').readlines()]


    stats={}
    ##get optimal number of edges
    opt_g=networkx.from_edgelist([(a[0],a[1]) for a in opt_tab])
    dum_g=networkx.from_edgelist([(a[0],a[1]) for a in dummy_tab])

    stats['numEdges'] = len(opt_tab)
    ##get number of trees (dummy - opt)
    stats['numTrees'] = networkx.number_connected_components(opt_g) #len(dummy_tab) - len(opt_tab)

    conn_comps=networkx.connected_components(opt_g)

    ##number of nodes
    nodes=opt_g.nodes()
    stats['numNodes']=len(nodes)

    ##get tree sizes
    stats['treeSizes']=','.join([str(len(c)) for c in sorted(conn_comps, key=len, reverse=True)])

    ##find ubiquitin?
    if 'UBC' in nodes:
    	stats['hasUBC']='True'
    else:
	stats['hasUBC']='False'

    return stats

##lots of for loops to query everything involved, and then read in file to
##collect stats:
##1- individual number of edges#
##2- number of trees
##3- whther or not ubiquitin was included

##goal it to write it all into tab-delimited file so we can
##analyze it in R.
ofile=open('forestStats.txt','w')
ofile.write('inputDataType\tpatientCoverage\tmu\tbeta\tw\tnumEdges\tnumTrees\tnumNodes\ttreeSizes\thasUBC\n')

for a in allfiles.keys():
    input_path = "../2016-03-18/"+a
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
                    opt_file=os.path.join(input_path,out_label+'_optimalForest.sif')
                    dum_file=os.path.join(input_path,out_label+'_dummyForest.sif')
                    if os.path.exists(opt_file) and os.path.exists(dum_file):
                        net_stats=netStats(opt_file,dum_file)
                        ofile.write("%s\t%s\t%f\t%d\t%f\t%d\t%d\t%d\t%s\t%s\n"%(a,prefix,mu,beta,w,net_stats['numEdges'],\
                                                                            net_stats['numTrees'],net_stats['numNodes'],\
                                                                            net_stats['treeSizes'],net_stats['hasUBC']))

ofile.close()

syn.store(synapseclient.File('forestStats.txt',parentId='syn5804586'),activityName='evaluatedForestAndParams',used = {'url':'https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-21/evalMoreForestNetworks.py'})
