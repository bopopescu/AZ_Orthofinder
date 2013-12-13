from os import chdir
import subprocess


def orthomcl_adjust_fasta(faa_dirpath):
    chdir(faa_dirpath)
    subprocess.call('orthomclAdjustFasta')


if __name__ == '__main__':
    orthomcl_adjust_fasta('../data/Escherichia_coli_K-12_proteomes')