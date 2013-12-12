from Bio import SeqIO
from cached_entrez import efetch_multiple

species_name = 'Yersinia pestis biovar Microtus complete'
genomes_dirpath = efetch_multiple(species_name)


