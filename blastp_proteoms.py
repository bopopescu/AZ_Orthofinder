from os import chdir
from os.path import join, realpath, dirname
from subprocess import call
from Bio.Blast.Applications import NcbiblastpCommandline
from fetch_proteomes import fetch_proteomes

file_dirpath = dirname(realpath(__file__))
orthomcl_dir = join(file_dirpath, 'orthomcl_software/bin')


def blastp_all_vs_all(fasta_dirpath):
    good_fasta_fpath = join(fasta_dirpath, 'goodProteins.fasta')
    chdir(fasta_dirpath)
    call([join(orthomcl_dir, 'orthomclFilterFasta'), fasta_dirpath, '10', '20'])
    chdir(file_dirpath)

    blastp_commandline = NcbiblastpCommandline(
        query=good_fasta_fpath,
        subject=good_fasta_fpath,
        outfmt=0,
        out=join(fasta_dirpath, 'blasted.xml'))
    print blastp_commandline
    stdout, stderr = blastp_commandline()


if __name__ == '__main__':
    species_name = 'Escherichia coli K-12'
    proteomes_dirpath = fetch_proteomes(species_name)
    blastp_all_vs_all(proteomes_dirpath)