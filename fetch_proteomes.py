from itertools import izip, count
from sys import stderr
from os import makedirs, chdir
from os.path import join, isdir
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from cached_entrez import efetch_multiple


proteomes_root_dirpath = '../data'


def fetch_proteomes(species_name):
    species_dirpath = join(
        proteomes_root_dirpath,
        species_name.replace(' ', '_'))
    proteomes_dirpath = join(species_dirpath, 'proteomes')

    if isdir(proteomes_dirpath):
        return proteomes_dirpath

    fasta_fpaths, gb_fpaths = efetch_multiple(species_name)

    i = 1
    for fasta_fpath, gb_fpath in izip(fasta_fpaths, gb_fpaths):
        rec = SeqIO.read(gb_fpath, 'genbank')
        features = [f for f in rec.features if f.type == 'CDS']
        print '%s: translating %d features' % (rec.id, len(features))

        genome_seq = SeqIO.read(fasta_fpath, 'fasta')
        if not isdir(proteomes_dirpath):
            makedirs(proteomes_dirpath)

        orthomcl_taxoncode = ''.join([w[0].lower() for w in species_name.split()[:2]]) + str(i)

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

            translation = Seq(qs.get('translation', [None])[0], generic_protein)
            if not translation:
                print >> stderr, '    Error: no translation for CDS'
            print '   ', translation

            # Translate ourselves
            #trans_table = int(qs.get('transl_table', [11])[0])
            #my_translation = f.extract(genome_seq).seq.translate(table=trans_table, to_stop=True, cds=True)
            #print my_translation
            #if translation:
            #    assert str(my_translation) == str(translation)

            proteins.append(SeqRecord(seq=translation, id=orthomcl_taxoncode + '|' + protein_id,
                                      description=protein_descripton))
            print

        if proteins:
            fpath = join(proteomes_dirpath, orthomcl_taxoncode + '.fasta')
            SeqIO.write(proteins, fpath, 'fasta')
            print 'Whitten to ' + fpath
            i += 1
        print

    return species_dirpath


if __name__ == '__main__':
    species_dirpath = fetch_proteomes('Escherichia coli K-12')