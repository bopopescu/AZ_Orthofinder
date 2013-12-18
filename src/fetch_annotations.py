from os import mkdir, listdir, remove
from os.path import join, isdir, isfile, splitext
from shutil import rmtree
from Bio import Entrez, Blast, SeqIO
Entrez.email = 'vladislav.sav@gmail.com'

import logging
log = logging.getLogger('orthofinder')

genbank_ext = 'gb'


#def efetch(db, id, rettype, **kwargs):
#    fpath = join(ref_dir, id + '.' + rettype)
#    if not isfile(fpath):
#        fetch_handle = Entrez.efetch(db, id=id, retmode='text', **kwargs)
#        record = Entrez.read(fetch_handle)
#        print 'efetch', record.id
#        with open(fpath, 'w') as f:
#            f.write(Entrez.read(fetch_handle))
#    return fpath


def fetch_annotations(save_dir, species_names, clip=None):
    if not species_names:
        log.error('No species names')
        return 1

    if not isdir(save_dir):
        mkdir(save_dir)

    for species_i, species_name in enumerate(species_names):
        species_name = species_name.strip()
        if species_name == '' or species_name[0] == '#':
            continue

        term = '%s[Organism] AND (complete genome[Title] OR complete sequence[Title])' % species_name
        log.info('   Query: %s' % term)
        search_handle = Entrez.esearch(db='nuccore', retmax=clip, term=term)
        ids = Entrez.read(search_handle)['IdList']

        if ids == []:
            log.info('No references are found.')
            ids = raw_input('Put reference ids manually:').split()
            if ids == []:
                log.error('No references :(')
                return 1

        log.info('   Found %s' % ids)

        for i, id in enumerate(ids):
            log.info('   Fetching %s...' % id)

            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype=genbank_ext)
            gb_fpath = join(save_dir, str(species_i) + '_' + str(i) + '_' + id + '.' + genbank_ext)
            with open(gb_fpath, 'w') as file:
                file.write(fetch_handle.read())

            rec = SeqIO.read(gb_fpath, genbank_ext)
            org_name = rec.annotations['organism']
            log.info('       ' + rec.annotations['organism'])
            if 'plasmid' in org_name:
                remove(gb_fpath)
            else:
                log.info('       saved %s' % gb_fpath)
            log.info('')

            #fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype=fasta_ext, **kwargs)
            #fasta_fpath = join(dirpath, str(i) + '_' + id + '.' + fasta_ext)
            #fasta_fpaths.append(fasta_fpath)
            #with open(fasta_fpath, 'w') as file:
            #    file.write(fetch_handle.read())
            #    print 'and .fasta'

    return 0


#def esearch(db, term, dirname, rettype):
#    dirpath = join(cache_dirpath, dirname)
#    fpaths = []
#
#    if not isdir(dirpath):
#        with Entrez.esearch(db, term, usehistory='y') as search_handle:
#            search_results = Entrez.read(search_handle)
#        print 'ids found:', search_results['IdList']
#
#        #batch_size = 10  # sequences to be downloaded at once
#        #for start in range(0, int(search_results['Count']), batch_size):
#        #    fetch_handle = Entrez.efetch(
#        #        db=db, rettype=rettype, retmode='text',
#        #        retstart=start, retmax=batch_size,
#        #        webenv=search_results['WebEnv'],
#        #        query_key=search_results['QueryKey'])
#
#        if (search_results['Count']) == 0
#        for id in search_results['IdList']:
#
#
#
#            data = fetch_handle.read()
#            fpaths = dict((id, join(dirpath, id + '.' + rettype)) for id in search_results['IdList'])
#
#        batch_size = 10
#        with open(
