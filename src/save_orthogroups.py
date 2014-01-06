from itertools import count, izip
from os import listdir
from os.path import join

from Bio import SeqIO

import logging
import config
log = logging.getLogger(config.log_fname)

#Gene = namedtuple('Gene', 'protein strain gi gene locus product description')

def save_orthogroups(annotations_dir, mcl_output, out, out_nice):
    mcl_f = open(mcl_output)
    out_f = open(out, 'w')
    nice_f = open(out_nice, 'w')

    strains = dict()
    max_lengths = count(0)
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

    for i, line in izip(count(1), mcl_f):
        print >> out_f, 'Orthogroup %d' % i
        print >> nice_f, 'Orthogroup %d' % i

        for gene in line.split():
            taxon_id, prot_id = gene.split('|')
            for l, val in zip(max_lengths, strains[taxon_id][prot_id][:-1]):
                print >> out_f, str(val) + '\t',
                print >> nice_f, str(val) + ' ' * (l - len(str(val))) + '\t',

            print >> out_f, str(strains[taxon_id][prot_id][-1])
            print >> nice_f, str(strains[taxon_id][prot_id][-1])

        print >> out_f
        print >> nice_f

    mcl_f.close()
    out_f.close()
    nice_f.close()

    return 0