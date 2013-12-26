from itertools import count, izip
from os import listdir
from os.path import join

from Bio import SeqIO

import logging
import config
log = logging.getLogger(config.log_fname)

#Gene = namedtuple('Gene', 'protein strain gi gene locus product description')

#TODO: make output columns align vertically, nicer!
def save_orthogroups_gff(annotations_dir, mcl_output, out):
    mcl_f = open(mcl_output)
    out_f = open(out, 'w')

    strains = dict()
    gb_files = [join(annotations_dir, fname)
                for fname in listdir(annotations_dir) if fname[0] != '.']
    for fname in gb_files:
        log.debug('   Reading ' + fname)

        rec = SeqIO.read(fname, 'genbank')
        strain = rec.annotations['source'] or rec.name
        gi = rec.annotations['gi'] or rec.id or 'NA'
        locus = rec.name
        description = rec.description

        genes_by_protid = dict()
        for feature in rec.features:
            if feature.type == 'CDS':
                qs = feature.qualifiers
                prot_id = qs.get('protein_id', ['NA'])[0]
                gene_id = qs.get('gene', ['NA'])[0]
                product = qs.get('product', ['NA'])[0]

                genes_by_protid[prot_id] = [prot_id, strain, gi, gene_id,
                                            locus, product, description]
        strains[rec.id] = genes_by_protid

    for i, line in izip(count(1), mcl_f):
        print >> out_f, 'Orthogroups %d' % i
        for gene in line.split():
            taxon_id, prot_id = gene.split('|')
            print >> out_f, '\t'.join(strains[taxon_id][prot_id])
        print >> out_f

    mcl_f.close()
    out_f.close()

    return 0