from os import mkdir, remove
from os.path import join, isdir, basename
from ftplib import FTP
from utils import read_list
from ftp_proxy import setup_http_proxy
from Bio import Entrez, SeqIO
Entrez.email = 'vladislav.sav@gmail.com'

import config
import logging
log = logging.getLogger(config.log_fname)


def fetch_annotations_for_species_from_ftp(save_dir, species_names, proxy=None, clip=None):
    if not species_names:
        log.error('   No species names')
        return 1

    if proxy:
        setup_http_proxy(*proxy.split(':'))

    if not isdir(save_dir):
        mkdir(save_dir)

    ftp = FTP('ftp.ncbi.nih.gov')
    log.debug('   Logging in ' + ftp.login())
    ftp.cwd('genomes/Bacteria')
    for i, dirname in enumerate(ftp.nlst()):
        for sp_i, sp in enumerate(species_names):
            sp = sp.replace(' ', '_').strip()
            if sp == '' or sp[0] == '#':
                continue

            if dirname.startswith(sp):
                log.debug('   Scanning ' + dirname)
                for remote_fpath in ftp.nlst(dirname):
                    if remote_fpath.endswith('.gbk'):
                        dest_fname = str(i) + '_' + str(sp_i) + '_' + basename(remote_fpath)
                        dest_fpath = join(save_dir, dest_fname)
                        ftp.retrbinary('RETR ' + remote_fpath, open(dest_fpath, 'wb').write)

                        rec = SeqIO.read(dest_fpath, 'gb')
                        #print '       Definition: ' + rec.description
                        if 'plasmid' in rec.description:
                            remove(dest_fpath)
                        else:
                            log.info('       saved ' + dest_fpath)
                            log.info('')
    return 0


def fetch_annotations_for_ids(annotations_dir, ref_ids):
    if not isdir(annotations_dir):
        mkdir(annotations_dir)

    if ref_ids == []:
        log.info('   No references have been found.')
        ids = raw_input('   Put reference ids manually:').split()
        if ids == []:
            log.error('   No references :(')
            return 1

    log.info('   IDs: %s' % ', '.join(ref_ids))

    for i, id in enumerate(ref_ids):
        log.info('   Fetching annotations for %s...' % id)

        try:
            fetch_handle = Entrez.efetch(db='nucleotide', id=id,
                                         retmode='text', rettype='gbwithparts')
            gb_fpath = join(annotations_dir, id + '.gb')
            with open(gb_fpath, 'w') as file:
                file.write(fetch_handle.read())

            rec = SeqIO.read(gb_fpath, 'genbank')
            genes_number = len([f for f in rec.features if f.type == 'CDS'])
            log.info('       ' + rec.description)
            log.info('       %d genes found.' % genes_number)
            log.info('       saved %s' % gb_fpath)
            log.info('')

        except KeyboardInterrupt, e:
            return 1

    return 0


def fetch_annotations_species_name_entrez(save_dir, species_names, clip=None):
    if not species_names:
        log.error('   No species names')
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
            log.info('   No references are found.')
            ids = raw_input('   Put reference ids manually:').split()
            if ids == []:
                log.error('   No references :(')
                return 1

        log.info('   IDs: %s' % ', '.join(ids))

        for i, id in enumerate(ids):
            log.info('   Fetching %s...' % id)

            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text',
                                         rettype='gbwithparts')
            gb_fpath = join(save_dir, str(species_i) + '_' + str(i) + '_' + id + '.' + genbank_ext)
            with open(gb_fpath, 'w') as file:
                file.write(fetch_handle.read())

            rec = SeqIO.read(gb_fpath, 'gb')
            log.info('       Organism: ' + rec.annotations['organism'])
            log.info('       Definition: ' + rec.description)
            if 'plasmid' in rec.description:
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


#def efetch(db, id, rettype, **kwargs):
#    fpath = join(ref_dir, id + '.' + rettype)
#    if not isfile(fpath):
#        fetch_handle = Entrez.efetch(db, id=id, retmode='text', **kwargs)
#        record = Entrez.read(fetch_handle)
#        print 'efetch', record.id
#        with open(fpath, 'w') as f:
#            f.write(Entrez.read(fetch_handle))
#    return fpath


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
