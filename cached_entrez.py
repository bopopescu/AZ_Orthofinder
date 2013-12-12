from os import mkdir, listdir
from os.path import join, isdir, isfile
from Bio import Entrez, Blast, SeqIO
Entrez.email = 'vladislav.sav@gmail.com'


cache_dirpath = 'cache'


def efetch(db, id, rettype, **kwargs):
    fpath = join(cache_dirpath, id + '.' + rettype)
    if not isfile(fpath):
        fetch_handle = Entrez.efetch(db, id=id, retmode='text', **kwargs)
        record = Entrez.read(fetch_handle)
        print 'efetch', record.id
        with open(fpath, 'w') as f:
            f.write(Entrez.read(fetch_handle))
    return fpath


def efetch_multiple(species_name, clip=5, **kwargs):
    dirpath = join(cache_dirpath, species_name.replace(' ', '_'))
    if not isdir(dirpath) or not listdir(dirpath):
        if not isdir(dirpath):
            mkdir(dirpath)

        search_handle = Entrez.esearch(db='nucleotide', term=species_name)
        ids = Entrez.read(search_handle)['IdList'][:clip]

        fasta_fpaths = dict((id, join(dirpath, id + '.fasta')) for id in ids)
        gb_fpaths = dict((id, join(dirpath, id + '.gb')) for id in ids)

        for id in ids:
            print '  fetching ' + id
            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype='gb', **kwargs)
            with open(gb_fpaths[id], 'w') as file:
                file.write(fetch_handle.read())

            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype='fasta', **kwargs)
            with open(fasta_fpaths[id], 'w') as file:
                file.write(fetch_handle.read())

    return dirpath




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
