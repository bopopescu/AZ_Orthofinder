from itertools import count, izip, repeat
from os import listdir, mkdir, rmdir
from os.path import join, isdir, splitext, basename
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


def save_orthogroups(assembly_proteins, annotations, mcl_output,
                     out, out_nice, out_short, assembly_singletones):
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


    # Assembly genes
    assembly_name = None
    assembly_proteins_recs = dict()
    if assembly_proteins:
        log.debug('   Reading ' + assembly_proteins)
        assembly_name = splitext(basename(assembly_proteins))[0]
        strain_id = assembly_name
        strain = assembly_name

        genes_by_protid = dict()
        for rec in SeqIO.parse(assembly_proteins, 'fasta'):
            prot_id = rec.id.split('|')[1]
            locus = 'NA'
            description = rec.description
            gene_id = 'NA'
            product = 'NA'
            assembly_proteins_recs[prot_id].append(rec)

            genes_by_protid[prot_id] = \
                [strain, prot_id, gene_id, locus, product, description]

            max_lengths = map(max, zip(max_lengths, map(len, genes_by_protid[prot_id][:-1])))

        strains[strain_id] = genes_by_protid


    # Other genes
    for fname in gb_files:
        log.debug('   Reading ' + fname)

        #if assembly_name and '.' in fname and splitext(basename(fname))[0] == assembly_name:
        #    strain_id = assembly_name
        #    strain = assembly_name
        #
        #    a = SeqIO.parse(fname, 'genbank')
        #    b = list(a)
        #
        #    for rec in SeqIO.parse(fname, 'genbank'):
        #        locus = rec.name
        #        description = rec.definition
        #
                #        genes_by_protid = dict()
        #
        #        for feature in rec.features:
        #            if feature.type == 'CDS':
        #                qs = feature.qualifiers
        #                prot_id = qs.get('protein_id', ['NA'])[0]
        #                gene_id = qs.get('gene', ['NA'])[0]
        #                product = qs.get('product', ['NA'])[0]
        #
        #                genes_by_protid[prot_id] = \
        #                    [strain, prot_id, gene_id, locus, product, description]
        #
        #                max_lengths = map(max, zip(max_lengths, map(len, genes_by_protid[prot_id][:-1])))
        #
        #        strains[strain_id] = genes_by_protid

        #else:
        try:
            rec = SeqIO.read(fname, 'genbank')
            strain_id = rec.id
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

        strains[strain_id] = genes_by_protid


    # Writing result files
    singletone_assembly_recs = []

    with open(mcl_output) as mcl_f, \
         open(out, 'w') as out_f, \
         open(out_nice, 'w') as nice_f, \
         open(assembly_singletones, 'w') as singletones_f:

        genes_number = 0
        i = 0
        for line in mcl_f:
            i += 1
            print >> out_f, 'Orthogroup %d' % i
            print >> nice_f, 'Orthogroup %d' % i

            known_genes_in_group = []
            for gene in line.split():
                genes_number += 1
                taxon_id, prot_id = gene.split('|')

                if taxon_id != assembly_name:
                    known_genes_in_group.append(prot_id)

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

            if assembly_name and known_genes_in_group == []:
                group = []
                for gene in line.split():
                    taxon_id, prot_id = gene.split('|')
                    group.append((i, assembly_proteins_recs[prot_id]))
                    #print >> singletones_f, '\t'.join(iter(strains[taxon_id][prot_id]))
                singletone_assembly_recs.append(group)

    log.info('   Saved %d groups, totally containing %d genes.' % (i, genes_number))

    for i, group in singletone_assembly_recs:
        SeqIO.write(group, splitext(assembly_singletones)[0] + '_group_' + str(i) + '.fasta', 'fasta')

    log.info('   Saved %d proteins, totally containing %d genes.' % (i, genes_number))

    return 0








