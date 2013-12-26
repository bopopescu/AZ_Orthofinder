from os.path import basename, join, relpath

from Workflow import Step
from save_orthogroups_gff import save_orthogroups_gff
from make_proteomes import make_proteomes
from fetch_annotations import fetch_annotations_for_species_from_ftp, fetch_annotations_for_ids
from config import orthomcl_config, orthomcl_bin_dir, BLAST_DBSIZE
from clean_db import clean_db

proteomes_dir       = 'proteomes'
annotations_dir     = 'annotations'
intermediate_dir    = 'intermediate'
sql_log             = 'intermediate/log.sql'
good_proteins       = 'intermediate/good_proteins.fasta'
poor_proteins       = 'intermediate/poor_proteins.fasta'
blast_db            = 'intermediate/blastdb'
blast_out           = 'intermediate/blasted.tsv'
similar_sequences   = 'intermediate/similar_sequences.txt'
pairs_log           = 'intermediate/orthomclpairs.log'
mcl_input           = 'intermediate/mcl_input'
mcl_output          = 'intermediate/mcl_output'
pairs_dir           = 'intermediate/pairs'
pairs_orthologs     = 'intermediate/pairs/potentialOrthologs.txt'
pairs_inparalogs    = 'intermediate/pairs/potentialInparalogs.txt'
pairs_coorthologs   = 'intermediate/pairs/potentialCoorthologs.txt'
groups_file         = 'groups.txt'
singletons_file     = 'singletons.txt'
orthogroups_file    = 'orthogroups.txt'


with open(orthomcl_config) as f:
    conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')
    ortholog_table = conf['orthologTable'].strip()
    in_paralog_table = conf['inParalogTable'].strip()
    coortholog_table = conf['coOrthologTable'].strip()
    similar_sequeces_table = conf['similarSequencesTable'].strip()
    inter_taxon_match_view = conf['interTaxonMatchView'].strip()
    best_hit_table = 'BestHit'
    best_hit_taxon_score_table = 'BestQueryTaxonScore'


def step_fetching_annotations_for_species(specied_list, proxy):
    return Step(
        'Fetching annotations',
         cmd=fetch_annotations_for_species_from_ftp,
         prod_files=[annotations_dir],
         parameters=[annotations_dir,
                     specied_list,
                     proxy])

def step_fetch_annotations_for_ids(ids_list):
    return Step(
        'Fetching annotations',
         cmd=fetch_annotations_for_ids,
         prod_files=[annotations_dir],
         parameters=[annotations_dir,
                     ids_list])

def step_make_proteomes(annotations_dir=annotations_dir):
    return Step(
        'Preparing proteomes',
         cmd=make_proteomes,
         req_files=[annotations_dir],
         prod_files=[proteomes_dir],
         parameters=[annotations_dir, proteomes_dir])

def filter_proteomes():
    return Step(
        'Filtering proteomes',
         cmd=join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
         req_files=[proteomes_dir],
         prod_files=[good_proteins,
                     poor_proteins],
         parameters=[proteomes_dir,
                     10, 20,
                     good_proteins,
                     poor_proteins])

def make_blast_db():
    return Step(
        'Making blast database',
         cmd='makeblastdb',
         req_files=[good_proteins],
         prod_files=[blast_db + '.' + ext for ext in ['phr', 'pin', 'psq']],
         parameters=[
            '-in', good_proteins,
            '-input_type', 'fasta',
            '-out', blast_db,
            '-dbtype', 'prot'],
         stdout='log')

def blast(threads):
    return Step(
        'Blasting',
         cmd='blastp',
         req_files=[good_proteins],
         prod_files=[blast_out],
         parameters=[
            '-query', good_proteins,
            '-db', blast_db,
            '-out', blast_out,
            '-outfmt', 6,  # tabular
            '-evalue', 1e-5,
            '-num_descriptions', 10000,  # don't care value
            '-num_alignments', 10000,  # don't care value
            '-num_threads', threads,
            '-dbsize', BLAST_DBSIZE])

def parse_blast_restults():
    return Step(
        'Parsing blast results',
         cmd=join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
         req_files=[proteomes_dir,
                    blast_out],
         prod_files=[similar_sequences],
         parameters=[blast_out,
                     proteomes_dir],
         stdout=similar_sequences)

def clean_database(suffix):
    return Step(
        'Cleaning database',
         cmd=clean_db,
         parameters=[suffix])

def install_schema(suffix):
    return Step(
        'Installing schema',
         cmd=join(orthomcl_bin_dir, 'orthomclInstallSchema.pl'),
         req_files=[orthomcl_config],
         prod_tables=[
             ortholog_table + suffix,
             in_paralog_table + suffix,
             coortholog_table + suffix,
             similar_sequeces_table + suffix,
             inter_taxon_match_view + suffix],
         parameters=[
             orthomcl_config,
             sql_log,
             suffix],
         stderr='log')

def load_blast_results(suffix):
    return Step(
        'Loading blast results into the database',
         cmd=join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
         req_files=[orthomcl_config,
                    similar_sequences],  # and initialized database
         prod_files=[],  # loads blast results into the db
         parameters=[orthomcl_config,
                     similar_sequences,
                     suffix],
         stderr='log')

def find_pairs(suffix):
    return Step(
        'Finding pairs',
         cmd=join(orthomcl_bin_dir, 'orthomclPairs.pl'),
         req_files=[orthomcl_config],
         req_tables=[in_paralog_table + suffix,
                     ortholog_table + suffix,
                     coortholog_table + suffix],
         prod_files=[],  # populates InParalog, Ortholog, CoOrtholog
         parameters=[orthomcl_config, pairs_log, 'cleanup=no', 'suffix=' + suffix],
         stderr='log')

def dump_pairs_to_files(suffix):
    return Step(
        'Dump pairs files',
         cmd=join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
         req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
         req_tables=[in_paralog_table + suffix,
                     ortholog_table + suffix,
                     coortholog_table + suffix],
         prod_files=[mcl_input,
                     pairs_dir,
                     pairs_orthologs,
                     pairs_inparalogs,
                     pairs_coorthologs],
         parameters=[orthomcl_config,
                     relpath(mcl_input, intermediate_dir),
                     intermediate_dir,
                     suffix],
         stderr='log')

def mcl(inflation=1.5):
    return Step(
        'MCL',
         cmd='mcl',
         req_files=[mcl_input],
         prod_files=[mcl_output],
         parameters=[mcl_input,
                     '--abc',
                     '-I', str(inflation),
                     '-o', mcl_output],
         stderr='log',
         stdout='log')

def step_save_orthogroups(annotations_dir=annotations_dir):
    return Step(
        'Saving orthogroups',
        cmd=save_orthogroups_gff,
        req_files=[mcl_output],
        prod_files=[orthogroups_file],
        parameters=[annotations_dir, mcl_output, orthogroups_file])

def groups_to_files(prefix, start_id):
    return Step(
        'MCL groups to files',
         cmd=join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
         req_files=[mcl_output],
         prod_files=[groups_file],
         parameters=[prefix + '_', start_id],
         stdin=mcl_output,
         stdout=groups_file)

def signletones_to_files():
    return Step(
        'MCL singletones to files',
         cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
         req_files=[good_proteins, groups_file],
         prod_files=[singletons_file],
         parameters=[good_proteins, groups_file],
         stdout=singletons_file)