from Bio import SeqIO
i = 0
for rec in SeqIO.parse('/Users/vladsaveliev/Dropbox/bio/az/EcoliKpneumo/intermediate/good_proteins.fasta', 'fasta'):
    i += 1

print i
