from itertools import count, izip
from os import listdir, mkdir
from os.path import join, isdir

from Bio import SeqIO, Entrez
from fetch_annotations import fetch_annotations_for_ids
Entrez.email = 'vladislav.sav@gmail.com'

import logging
import config
log = logging.getLogger(config.log_fname)

#Gene = namedtuple('Gene', 'protein strain gi gene locus product description')

def save_orthogroups(annotations_dir, mcl_output, out, out_nice):
    strains = dict()
    max_lengths = count(0)

    if not isdir(annotations_dir) or not listdir(annotations_dir):
        log.info('')
        log.info('   Fetching annotations from Genbank.')

        ids = set()
        with open(mcl_output) as mcl_f:
            for line in mcl_f:
                for gene in line.split():
                    taxon_id, _ = gene.split('|')
                    ids.add(taxon_id)

        fetch_annotations_for_ids(annotations_dir, ids)

    gb_files = [join(annotations_dir, fname)
                for fname in listdir(annotations_dir) if fname[0] != '.']

    for fname in gb_files:
        log.debug('   Reading ' + fname)

        rec = SeqIO.read(fname, 'genbank')
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
                for i, line in izip(count(1), mcl_f):
                    print >> out_f, 'Orthogroup %d' % i
                    print >> nice_f, 'Orthogroup %d' % i

                    for gene in line.split():
                        taxon_id, prot_id = gene.split('|')
                        if taxon_id not in strains:
                            log.error('   No annotations for ' + taxon_id + ' in ' + annotations_dir)
                            return 1

                        for l, val in zip(max_lengths, strains[taxon_id][prot_id][:-1]):
                            print >> out_f, str(val) + '\t',
                            print >> nice_f, str(val) + ' ' * (l - len(str(val))) + '\t',

                        print >> out_f, str(strains[taxon_id][prot_id][-1])
                        print >> nice_f, str(strains[taxon_id][prot_id][-1])

                    print >> out_f
                    print >> nice_f

    return 0