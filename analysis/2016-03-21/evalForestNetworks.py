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
beta_range = [1,5,10]
w_range = [1,2,3]

edge_file = os.path.join(omics_int_dir,"data/iref_mitab_miscore_2013_08_12_interactome.txt")
conf_file = "conf.txt"
D = 5

'''
Get network statistics for a particular .sif file to try to determine best
one to use!
'''
def netStats(opt_fname,dummy_fname):
    opt.tab=[a.split('\tpp\t') for a in open(opt_fname,'r').readLines()]
    dummy.tab=[a.split('\tpp\t') for a in open(dummy_fname,'r').readLines()]

    stats={}
    ##get optimal number of edges
    opt.g=networkx.from_edgelist(opt.tab)
    dum.g=networkx.from_edgelist(dumm.tab)

    stats['numEdges'] = len(opt.tab)
    ##get number of trees (dummy - opt)
    stats['numTrees'] = len(dummy.tab) - len(opt.tab)

    conn_comps=networkx.connected_components(opt.g)
    print stats['numTrees'],conn_comps

    ##number of nodes
    nodes=set()
    [nodes.add([a[0],a[1]]) for a in opt.tab]
    stats['numNodes']=len(nodes)

    ##get tree sizes
    stats['treeSize']=','.join([len(c) for c in sorted(conn_comps, key=len, reverse=True)])

    ##find ubiquitin?
    stats['hasUBC']='UBC'%in%nodes

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
                    net_stats=netStats(opt_file,dum_file)
                    ofile.write("%s\t%s\t%f\t%d\t%f\t%d\t%d\t%d\t%s\t%s\n"%(a,prefix,mu,beta,w,net_stats['numEdges'],net_stats['numTrees'],net_stats['numNodes'],net_stats['treeSizes'],net_stats['hasUBC'])
of.close()
