from os.path import join
from Bio import Entrez, SeqIO
Entrez.email = 'vladislav.sav@gmail.com'

import logging
log = logging.getLogger('orthofinder')

handle = Entrez.esearch(db="genome", term="Klebsiella pneumoniae")
record = Entrez.read(handle)
print record.keys()
if record['IdList']:
    handle = Entrez.efetch(db='genome', id=record['IdList'])
    rec = Entrez.read(handle)
    print rec
    print rec.keys()

#for id in record['IdList']:
#    print(id)
#    handle = Entrez.efetch(db="genome", id=id)
#    print(handle.read())
#
#for i, id in enumerate(record['IdList']):
#    print '   Fetching %s...' % id
#
#    fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype='gb')
#    gb_fpath = join(str(i) + '_test.gb')
#    with open(gb_fpath, 'w') as file:
#        file.write(fetch_handle.read())
#
#    rec = SeqIO.read(gb_fpath, 'gb')
#    org_name = rec.annotations['organism']
#    definition = rec.description
#    log.info('       Id: ' + rec.id)
#    log.info('       Organism: ' + org_name)
#    log.info('       Definition: ' + definition)
#    log.info('')
