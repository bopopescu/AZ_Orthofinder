from Bio import Entrez
Entrez.email = 'vladislav.sav@gmail.com'

handle = Entrez.esearch(db="pubmed", term="quast")
record = Entrez.read(handle)
for id in record['IdList']:
    handle = Entrez.efetch(db="pubmed", id=id)
    print(handle.read())





