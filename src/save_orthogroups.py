from genericpath import exists
from itertools import count, izip, repeat, chain
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


def get_assembly_genes(assembly_proteins_fpath, max_lengths):
    assembly_proteins_recs = dict()
    log.debug('   Reading ' + assembly_proteins_fpath)
    assembly_name = splitext(basename(assembly_proteins_fpath))[0]
    strain = assembly_name

    genes_by_protid = dict()
    for rec in SeqIO.parse(assembly_proteins_fpath, 'fasta'):
        prot_id = rec.id.split('|')[1]
        locus_tag = 'NA'
        description = rec.description
        gene_id = 'NA'
        product = 'NA'
        strain_id = assembly_name
        assembly_proteins_recs[rec.id] = rec

        genes_by_protid[prot_id] = \
            [strain, strain_id, prot_id, locus_tag, gene_id, product, description]

        max_lengths = map(max, zip(max_lengths, map(len, genes_by_protid[prot_id][:-1])))

    return genes_by_protid, assembly_proteins_recs, max_lengths


def get_reference_genes(fname, max_lengths):
    log.debug('   Reading ' + fname)
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
    description = rec.description

    genes_by_protid = dict()

    for feature in rec.features:
        if feature.type == 'CDS':
            qs = feature.qualifiers
            prot_id = qs.get('protein_id', ['NA'])[0]
            gene_id = qs.get('gene', ['NA'])[0]
            product = qs.get('product', ['NA'])[0]
            locus_tag = qs.get('locus_tag', ['NA'])[0]

            genes_by_protid[prot_id] = \
                [strain, strain_id, prot_id, locus_tag, gene_id, product, description]

            max_lengths = map(max, zip(max_lengths, map(len, genes_by_protid[prot_id][:-1])))

    return strain_id, genes_by_protid, max_lengths


def save_orthogroups(new_proteomes, annotations, mcl_output,
                     out, out_nice, out_short, assembly_singletones, singletone_dir):
    strains = dict()
    max_lengths = repeat(0)

    gb_files = []
    if isinstance(annotations, (list, tuple)):
        gb_files = annotations
    else:
        if isdir(annotations) and listdir(annotations):
            gb_files = [
                join(annotations, fname)
                for fname in listdir(annotations) if fname[0] != '.']
        #    if not isdir(annotations): mkdir(annotations)
        #    if __download(annotations, mcl_output) != 0:
        #        return 1

    if not gb_files:
        return save_compact(mcl_output, out)

    assembly_names = []
    new_proteins_recs = dict()

    for new_prot in new_proteomes:
        assembly_name = splitext(basename(new_prot))[0]
        assembly_names.append(assembly_name)

        genes, assembly_recs, max_lengths = \
            get_assembly_genes(new_prot, max_lengths)

        new_proteins_recs.update(assembly_recs)
        strains[assembly_name] = genes

    if new_proteins_recs:
        if exists(singletone_dir):
            rmtree(singletone_dir)
        if not isdir(singletone_dir):
            mkdir(singletone_dir)

    for fname in gb_files:
        strain_id, genes, max_lengths = get_reference_genes(fname, max_lengths)
        strains[strain_id] = genes

    with open(mcl_output) as mcl_f:
        groups_total = sum(1 for _ in mcl_f)

    with open(mcl_output) as mcl_f, \
         open(out, 'w') as out_f, \
         open(out_nice, 'w') as nice_f, \
         open(assembly_singletones, 'a') as singletones_f:

        singletone_assembly_recs = []

        gene_number = 0
        group_nunber = 0
        singletone_gene_number = 0
        singletone_group_number = 0

        for group_line in mcl_f:
            group_nunber += 1

            known_genes_in_this_group = []
            for rec_id in group_line.split():
                gene_number += 1
                taxon_id, prot_id = rec_id.split('|')

                if assembly_names and taxon_id not in assembly_names:
                    known_genes_in_this_group.append(prot_id)

                if taxon_id not in strains:
                    log.warn('   No annotations for "' + taxon_id + '"')
                    return 1

                for l, val in izip(chain([len(str(groups_total))], max_lengths),
                                   chain([group_nunber], strains[taxon_id][prot_id])):
                    out_f.write(str(val) + '\t')
                    nice_f.write(str(val) + ' ' * (l - len(str(val))) + '\t')
                out_f.write('\n')
                nice_f.write('\n')
            nice_f.write('\n')

            if new_proteins_recs and known_genes_in_this_group == []:
                group = []
                singletone_group_number += 1

                singletones_f.write('Group %d\n' % group_nunber)
                for rec_id in group_line.split():
                    singletone_gene_number += 1
                    group.append(new_proteins_recs[rec_id])
                    singletones_f.write(rec_id + ' ')
                singletones_f.write('\n')

                singletone_assembly_recs.append(group)

                a_singletone_filepath = join(
                    singletone_dir,
                    splitext(assembly_singletones)[0] + '_group_'
                    + str(group_nunber) + '.fasta')
                SeqIO.write(group, a_singletone_filepath, 'fasta')

    log.info('   Saved %d groups, totally containing %d genes.' % (group_nunber, gene_number))

    if singletone_assembly_recs:
        log.info('   Saved %d singletone groups for the assembly, totally containing %d genes.' %
                 (singletone_group_number, singletone_gene_number))

    return 0








