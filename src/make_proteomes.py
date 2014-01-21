from genericpath import isfile
from itertools import izip, count
from shutil import rmtree, copy, copyfile
import subprocess
from sys import stderr
from os import makedirs, chdir, mkdir, listdir, remove
from os.path import join, isdir, basename, splitext
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

import config
import logging
log = logging.getLogger(config.log_fname)


def adjust_proteomes(proteomes, proteomes_dir, prot_id_field):
    if not isdir(proteomes_dir):
        mkdir(proteomes_dir)

    for proteome in proteomes:
        prot_ids = set()
        records = []
        taxon_id, ext = splitext(basename(proteome))
        ext = '.fasta'
        for seq in SeqIO.parse(proteome, 'fasta'):
            fields = seq.id.replace('|', ' ').split()
            prot_id = fields[prot_id_field]
            if len(prot_id) > 30:
                prot_id = prot_id[:17] + '...' + prot_id[-10:]
            if prot_id in prot_ids:
                log.error('Fasta %s contains duplicate id: %s' % (proteome, prot_id))
                return 1
            prot_ids.add(prot_id)
            seq.id = taxon_id + '|' + prot_id
            records.append(seq)
        out = join(proteomes_dir, taxon_id + ext)
        if isfile(out):
            remove(out)
        SeqIO.write(records, out, 'fasta')

    return 0


def make_proteomes(annotations, proteomes_dir):
    gb_files = None
    if isinstance(annotations, (list, tuple)):
        gb_files = annotations

    elif isdir(annotations) and listdir(annotations):
        annotations_dir = annotations
        gb_files = [
            join(annotations_dir, fname)
            for fname in listdir(annotations_dir)
            if fname and fname[0] != '.']

    if not gb_files: log.error('   No references provided.')

    if not isdir(proteomes_dir): makedirs(proteomes_dir)

    #words = ' '.join(species_names).split()
    #speciescode = workflow_id[:4]  # ''.join(w[0].lower() for w in words[:3 - int(math.log10(ref_num))])

    i = 1
    for gb_fpath in gb_files:
        try:
            rec = SeqIO.read(gb_fpath, 'genbank')
        except:
            log.warning('   Cannot read proteins from ' + gb_fpath)
            continue

        features = [f for f in rec.features if f.type == 'CDS']
        log.info('   %s: translating %d features' % (rec.id, len(features)))

        taxoncode = rec.id

        proteins = []
        for f in features:
            qs = f.qualifiers
            protein_id = qs.get('protein_id', [None])[0]
            gene_id = qs.get('gene', [None])[0]
            if not protein_id:
                log.warn('   Warning: no protein_id for CDS')
                continue
            #if not gene_id:
            #    log.warn('   Warning: no gene_id for CDS')
            #    continue

            protein_descripton = rec.id + ' ' + rec.description + \
                                 ' ' + 'Gene ' + (gene_id or '<unknown>') + \
                                 '. Protein ' + protein_id
            #log.debug('   ' + protein_descripton)
            translation = None
            if qs.get('translation', [None]) is not None:
                translation = Seq(qs.get('translation', [None])[0], generic_protein)

            if not translation:
                # Translate ourselves
                # TODO: Fetch reference
                # read genome_seq
                #trans_table = int(qs.get('transl_table', [11])[0])
                #my_translation = f.extract(genome_seq).seq.translate(table=trans_table, to_stop=True, cds=True)
                #print my_translation
                #if translation:
                #    assert str(my_translation) == str(translation)

                # TODO: OR BETTER FETCH PROTEIN
                fetch_handle = Entrez.efetch(db='protein', id=protein_id,
                                             retmode='text', rettype='fasta')
                rec = SeqIO.read(fetch_handle, 'fasta')
                translation = rec.seq

            proteins.append(SeqRecord(seq=translation, id=taxoncode + '|' + protein_id,
                                      description=protein_descripton))
            #log.debug('')

        if proteins:
            fpath = join(proteomes_dir, taxoncode + '.fasta')
            SeqIO.write(proteins, fpath, 'fasta')
            log.info('   Written to ' + fpath)
            i += 1
        log.info('')

    return 0


#if __name__ == '__main__':
    #species_dirpath = make_proteomes('Escherichia coli K-12')