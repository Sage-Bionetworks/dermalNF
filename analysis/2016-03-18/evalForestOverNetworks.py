'''
Goal is to download files from synapse and iterate over different parameters
'''


import os,re
import synapseclient

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

##lots of for loops to query everything involved!
for a in allfiles.keys():
    output_path = "./"+a
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
                    conf_prep(mu,beta,D,w)
                    out_label = "%s_w%f_beta%d_D%d_mu%f" %(prefix,w,beta,D,mu)
                    cmd="python %s --prize %s --edge %s --conf conf.txt --msgpath ~/msgsteiner-1.1/msgsteiner --outpath %s --outlabel %s" %(os.path.join(omics_int_dir,'scripts/forest.py'),prize_file,edge_file,output_path,out_label)
                    print cmd
                    os.system(cmd)
