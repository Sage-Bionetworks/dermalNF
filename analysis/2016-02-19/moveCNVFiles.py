'''
Basic script to move around some files in the  CNV analysis directory
'''


import synapseclient
syn = synapseclient.login()

#first get all files in the CNV directory (that are not directories themselves)
sq=syn.query("select name,id from entity where parentId=='syn5049702'")

eids=[a['entity.id'] for a in sq['results'] if ('png' in a['entity.name'] or
                                                'pdf' in a['entity.name'] or
                                                'txt' in a['entity.name'])]

print eids

#all left should be gene/region analysis
for e in eids:
    sf=syn.get(e,downloadFile=False)
    print sf.name
    sf['parentId']='syn5669822'
    syn.store(sf)
