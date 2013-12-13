from sys import stderr
from os import makedirs
from os.path import join, isdir
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from cached_entrez import efetch_multiple


proteomes_root_dirpath = '../data'


def fetch_proteomes(species_name):
    fasta_fpaths, gb_fpaths = efetch_multiple(species_name)

    for fasta_fpath, gb_fpath in zip(fasta_fpaths, gb_fpaths):
        rec = SeqIO.read(gb_fpath, 'genbank')
        features = [f for f in rec.features if f.type == 'CDS']
        print '%s: translating %d features' % (rec.id, len(features))

        genome_seq = SeqIO.read(fasta_fpath, 'fasta')
        proteomes_dirpath = join(proteomes_root_dirpath,
                                 species_name.replace(' ', '_') + '_proteomes')
        if not isdir(proteomes_dirpath):
            makedirs(proteomes_dirpath)

        proteins = []
        for f in features:
            qs = f.qualifiers
            protein_id = qs.get('protein_id', [None])[0]
            gene_id = qs.get('gene', [None])[0]
            if not protein_id:
                print >> stderr, '    Warning: no protein_id for CDS'
                continue
            if not gene_id:
                print >> stderr, '    Warning: no gene_id for CDS'
                continue

            protein_descripton = '    ' + rec.id + ' ' + rec.description + \
                                 ' ' + 'Gene ' + gene_id + '. Protein ' + protein_id
            print protein_descripton

            annotation_translation = Seq(qs.get('translation', [None])[0], generic_protein)
            if not annotation_translation:
                print >> stderr, '    Error: no translation for CDS'

            print '   ', annotation_translation

            # Translate ourselves
            #trans_table = int(qs.get('transl_table', [11])[0])
            #translation = f.extract(genome_seq).seq.translate(table=trans_table, to_stop=True, cds=True)
            #print translation
            #if annotation_translation:
            #    assert str(translation) == str(annotation_translation)

            proteins.append(SeqRecord(seq=annotation_translation, id=protein_id, description=protein_descripton))
            print

        if proteins:
            fpath = join(proteomes_dirpath, rec.id + '.faa')
            SeqIO.write(proteins, fpath, 'fasta')
            print 'Whitten to ' + fpath
        print

if __name__ == '__main__':
    fetch_proteomes('Escherichia coli K-12')