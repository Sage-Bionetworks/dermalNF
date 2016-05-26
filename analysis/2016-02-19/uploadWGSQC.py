'''
Get WGS qc Analysis from 2016-02-09 and upload it to synapse
'''

import synapseclient
syn = synapseclient.login()


tabfile='../2016-02-09/allSamples.tab'
pyfile='../2016-02-09/allSamples.py'
pfile='../2016-02-09/allSamples.png'

parid='syn5669984'

for file in [tabfile,pyfile,pfile]:
    syn.store(synapseclient.File(file,parent=parid),activity=synapseclient.Activity('bcftools gtcheck',used='syn5555584'))
