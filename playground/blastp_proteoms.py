from os import chdir
from os.path import join, realpath, dirname
from subprocess import call
from Bio.Blast.Applications import NcbiblastpCommandline
from fetch_proteomes import fetch_proteomes
file_dirpath = dirname(realpath(__file__))


orthomcl_bin_dirpath = join(file_dirpath, 'orthomcl_software/bin')
config_fpath = join(file_dirpath, 'orthomcl.config')



def blastp_all_vs_all(species_dirpath):
    chdir(species_dirpath)

    good_fasta_fpath = 'goodProteins.fasta'

    blastp_commandline = NcbiblastpCommandline(
        query=good_fasta_fpath,
        subject=good_fasta_fpath,
        outfmt=6,  # tabular
        evalue=1e-5,
        #seg='"m S"',  # filter with seg
        num_descriptions=10000,  # don't care value
        num_alignments=10000,  # don't care value
        out='blasted_tabular')
    print blastp_commandline
    stdout, stderr = blastp_commandline()


def orthomcl(species_dirpath):
    for line in [l.strip() for l in '''
    orthomclFilterFasta proteomes_dir 10 20
    blastp -out blasted_tabular -outfmt 6 -num_descriptions 10000 -num_alignments 10000 -query goodProteins.fasta -evalue 1e-05 -subject goodProteins.fasta
    orthomclBlastParser blasted_tabular proteomes_dir >> similarSequences.txt
    orthomclLoadBlast orthomcl.config similarSequences.txt
    orthomclPairs orthomcl.config orthomclpairs.log cleanup=no
    orthomclDumpPairsFiles orthomcl.config
    mcl_software mclInput --abc -I 1.5 -o mclOutput
    orthomclMclToGroups az 1000 < mclOutput > groups.txt
    orthomclSingletons goodProteins.fasta groups.txt >> singletons.txt
    '''.split('\n') if l]:
        print line
        call(line.split())


if __name__ == '__main__':
    #species_name = 'Escherichia coli K-12'
    #species_dirpath = fetch_proteomes(species_name)
    #blastp_all_vs_all(species_dirpath)

    blastp_all_vs_all(join(file_dirpath, '../data/Escherichia_coli_K-12_2'))