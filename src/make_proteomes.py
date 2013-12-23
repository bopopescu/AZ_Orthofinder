from itertools import izip, count
from shutil import rmtree
from sys import stderr
from os import makedirs, chdir, mkdir, listdir
from os.path import join, isdir
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
import math

import config
import logging
log = logging.getLogger(config.log_fname)


def make_proteomes(gbk_dir, workflow_id, out_dir):
    gbk_files = [join(gbk_dir, fname) for fname in listdir(gbk_dir) if fname[0] != '.']
    ref_num = len(gbk_files)
    if ref_num == 0:
        log.error('   No references in ' + gbk_dir)
    if ref_num >= 1000:
        log.error('   Maximum 999 references are supported by OrthoMCL')

    if not isdir(out_dir):
        makedirs(out_dir)

    #words = ' '.join(species_names).split()
    speciescode = workflow_id[:4]  # ''.join(w[0].lower() for w in words[:3 - int(math.log10(ref_num))])

    for i, gb_fpath in izip(count(1), gbk_files):
        try:
            rec = SeqIO.read(gb_fpath, 'genbank')
        except:
            log.warning('   Cannot read proteins from ' + gb_fpath)
            continue

        features = [f for f in rec.features if f.type == 'CDS']
        log.info('   %s: translating %d features' % (rec.id, len(features)))

        taxoncode = speciescode + str(i)

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
            fpath = join(out_dir, taxoncode + '.fasta')
            SeqIO.write(proteins, fpath, 'fasta')
            log.info('   Written to ' + fpath)
        log.info('')

    return 0


#if __name__ == '__main__':
    #species_dirpath = make_proteomes('Escherichia coli K-12')