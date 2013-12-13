from os import mkdir, listdir
from os.path import join, isdir, isfile, splitext
from Bio import Entrez, Blast, SeqIO
Entrez.email = 'vladislav.sav@gmail.com'


cache_dirpath = 'cache'
fasta_ext = 'fasta'
genbank_ext = 'gb'


def efetch(db, id, rettype, **kwargs):
    fpath = join(cache_dirpath, id + '.' + rettype)
    if not isfile(fpath):
        fetch_handle = Entrez.efetch(db, id=id, retmode='text', **kwargs)
        record = Entrez.read(fetch_handle)
        print 'efetch', record.id
        with open(fpath, 'w') as f:
            f.write(Entrez.read(fetch_handle))
    return fpath


def efetch_multiple(species_name, clip=None, **kwargs):
    dirpath = join(cache_dirpath, species_name.replace(' ', '_'))

    if isdir(dirpath):
        return ([join(dirpath, f) for f in listdir(dirpath) if splitext(f)[1] == '.' + fasta_ext],
                [join(dirpath, f) for f in listdir(dirpath) if splitext(f)[1] == '.' + genbank_ext])

    else:
        mkdir(dirpath)
        fasta_fpaths, gb_fpaths = [], []

        term = '%s[Organism] AND (complete genome[Title] OR complete sequence[Title])' % species_name
        print 'Query: %s' % term
        search_handle = Entrez.esearch(db='nuccore', retmax=clip, term=term)
        ids = Entrez.read(search_handle)['IdList']
        print 'Found %s' % ids

        for i, id in enumerate(ids):
            print 'fetching %s...' % id,

            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype=genbank_ext, **kwargs)

            gb_fpath = join(dirpath, str(i) + '_' + id + '.' + genbank_ext)
            gb_fpaths.append(gb_fpath)
            with open(gb_fpath, 'w') as file:
                file.write(fetch_handle.read())
                print 'saved %s' % gb_fpath,

            fetch_handle = Entrez.efetch(db='nuccore', id=id, retmode='text', rettype=fasta_ext, **kwargs)

            fasta_fpath = join(dirpath, str(i) + '_' + id + '.' + fasta_ext)
            fasta_fpaths.append(fasta_fpath)
            with open(fasta_fpath, 'w') as file:
                file.write(fetch_handle.read())
                print ', saved %s' % fasta_fpath,

        return fasta_fpaths, gb_fpaths




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
