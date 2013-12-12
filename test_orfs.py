import os
from Bio import Entrez, Blast, SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from CachedEntrez import efetch, efetch_multiple

Entrez.email = 'vladislav.sav@gmail.com'

#handle = Entrez.esearch(db="Taxonomy", term="Homo sapiens")
#record = Entrez.read(handle)
#print(record['IdList'])
#
#id = record['IdList'][0]
#handle = Entrez.efetch(db='Taxonomy', id=id, retmode='xml')
#records = Entrez.read(handle)
#print(records[0]["Lineage"])

#gb_fname = 'data/NC_005816.gb'
#fasta_fname = 'data/NC_005816.fasta'


#records = SeqIO.read(fpath, 'gb')


#for strand, nuc in [(+1, record.seq),
#                    (-1, record.seq.reverse_complement())]:
#    for frame in [0, 1, 2]:
#        length = 3 * ((len(record) - frame) / 3)
#        for prot in nuc[frame : frame + length].translate(table=11).split('*'):
#            if len(prot) >= 100:
#                print("%s...%s - length %i, strand %i, frame %i" %
#                     (prot[:30], prot[-3:], len(prot), strand, frame))


#record = SeqIO.read(open(gb_fname), 'genbank')
#
#gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
#gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
#gd_feature_set = gd_track_for_features.new_set()
#
#for feature in record.features:
#    if feature.type != 'gene':
#        continue
#    if len(gd_feature_set) % 2 == 0:
#        color = colors.blue
#    else:
#        color = colors.lightblue
#    gd_feature_set.add_feature(feature, color=color, label=True)
#
#gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
#                fragments=4, start=0, end=len(record))
#
#gd_diagram.write("plasmid_linear.png", "PNG")
