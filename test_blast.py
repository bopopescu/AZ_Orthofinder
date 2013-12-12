import os
from Bio.Blast import NCBIWWW
from Bio import SeqIO

record = next(SeqIO.parse('datasets/ecoli.fasta', format='fasta'))

blast_result_fname = 'ecoli_blast.xml'

if not os.path.isfile(blast_result_fname):
    with NCBIWWW.qblast('blastn', 'nt', record.format('fasta')) as result_handle:
        with open(blast_result_fname, 'w') as f:
            f.write(result_handle.read())

#with open(blast_result_fname) as result_handle:

