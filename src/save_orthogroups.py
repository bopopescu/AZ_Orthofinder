from itertools import count, izip, repeat
from os import listdir, mkdir, rmdir
from os.path import join, isdir, splitext
from shutil import rmtree

from Bio import SeqIO, Entrez
from fetch_annotations import fetch_annotations_for_ids
Entrez.email = 'vladislav.sav@gmail.com'

import logging
import config
log = logging.getLogger(config.log_fname)

#Gene = namedtuple('Gene', 'protein strain gi gene locus product description')

#def __download(annotations, mcl_output):
#    log.info('   Fetching annotations from Genbank.')
#
#    ids = set()
#    with open(mcl_output) as mcl_f:
#        for line in mcl_f:
#            for gene in line.split():
#                taxon_id, prot_id = gene.split('|')
#                #h = Entrez.efetch(db='protein', id=id,
#                #                  retmode='text', rettype='gbwithparts')
#                ids.add(taxon_id)
#
#    return fetch_annotations_for_ids(annotations, ids)


def save_compact(mcl_output, out):
    with open(mcl_output) as mcl_f:
        with open(out, 'w') as out_f:
            for i, line in enumerate(mcl_f):
                out_f.write(str(i) + ' ' + line + '\n')
    return 0


def save_orthogroups(annotations, mcl_output, out, out_nice, out_short):
    strains = dict()
    max_lengths = count(0)

    gb_files = []
    if isinstance(annotations, (list, tuple)):
        gb_files = annotations
    else:
        if isdir(annotations) and listdir(annotations):
            gb_files = [join(annotations, fname)
                        for fname in listdir(annotations) if fname[0] != '.']
        #    if not isdir(annotations): mkdir(annotations)
        #    if __download(annotations, mcl_output) != 0:
        #        return 1

    if not gb_files:
        return save_compact(mcl_output, out)

    for fname in gb_files:
        log.debug('   Reading ' + fname)

        #if splitext(fname)
        try:
            rec = SeqIO.read(fname, 'genbank')
        except ValueError:
            log.error('   Could not read annotations from ' + fname)
            return 1
            #if isdir(annotations):
            #    rmtree(annotations)
            #    mkdir(annotations)
            #if __download(annotations, mcl_output) != 0:
            #    return 1
            #rec = SeqIO.read(fname, 'genbank')

        strain = rec.annotations['source'] or rec.name
        #gi = rec.annotations['gi'] or rec.id or 'NA'
        locus = rec.name
        description = rec.description

        genes_by_protid = dict()

        for feature in rec.features:
            if feature.type == 'CDS':
                qs = feature.qualifiers
                prot_id = qs.get('protein_id', ['NA'])[0]
                gene_id = qs.get('gene', ['NA'])[0]
                product = qs.get('product', ['NA'])[0]

                genes_by_protid[prot_id] = \
                    [strain, prot_id, gene_id, locus, product, description]

                max_lengths = map(max, zip(max_lengths, map(len, genes_by_protid[prot_id][:-1])))

        strains[rec.id] = genes_by_protid

    with open(mcl_output) as mcl_f:
        with open(out, 'w') as out_f:
            with open(out_nice, 'w') as nice_f:
                genes_number = 0
                i = 0
                for line in mcl_f:
                    i += 1
                    print >> out_f, 'Orthogroup %d' % i
                    print >> nice_f, 'Orthogroup %d' % i

                    for gene in line.split():
                        genes_number += 1
                        taxon_id, prot_id = gene.split('|')
                        if taxon_id not in strains:
                            log.warn('   No annotations for "' + taxon_id + '"')
                            return 1
                            #vals = repeat('NA')
                        else:
                            vals = iter(strains[taxon_id][prot_id])

                        for l, val in izip(max_lengths, vals):
                            print >> out_f, str(val) + '\t',
                            print >> nice_f, str(val) + ' ' * (l - len(str(val))) + '\t',

                        val = next(vals)
                        print >> out_f, val
                        print >> nice_f, val

                    print >> out_f
                    print >> nice_f

                log.info('   Saved %d groups, totally contating %d genes.' % (i, genes_number))

    return 0