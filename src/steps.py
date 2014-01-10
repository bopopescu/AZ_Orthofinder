from os.path import basename, join, relpath

from Workflow import Step, cmdline
from save_orthogroups import save_orthogroups
from make_proteomes import make_proteomes, adjust_proteomes
from fetch_annotations import fetch_annotations_for_species_from_ftp, fetch_annotations_for_ids
from config import orthomcl_config, orthomcl_bin_dir, BLAST_DBSIZE
from clean_db import clean_db

import config
import logging
log = logging.getLogger(config.log_fname)

proteomes_dir         = 'proteomes'
annotations_dir       = 'annotations'
intermediate_dir      = 'intermediate'
sql_log               = 'intermediate/log.sql'
predicted_proteins    = 'intermediate/predicted_proteins.fasta'
good_proteins         = 'intermediate/good_proteins.fasta'
poor_proteins         = 'intermediate/poor_proteins.fasta'
blast_db              = 'intermediate/blastdb'
blast_out             = 'intermediate/blasted.tsv'
similar_sequences     = 'intermediate/similar_sequences.txt'
pairs_log             = 'intermediate/orthomclpairs.log'
mcl_input             = 'intermediate/mcl_input'
mcl_output            = 'intermediate/mcl_output'
pairs_dir             = 'intermediate/pairs'
pairs_orthologs       = 'intermediate/pairs/potentialOrthologs.txt'
pairs_inparalogs      = 'intermediate/pairs/potentialInparalogs.txt'
pairs_coorthologs     = 'intermediate/pairs/potentialCoorthologs.txt'
groups_file           = 'groups.txt'
singletons_file       = 'singletons.txt'
orthogroups_file      = 'orthogroups.tsv'
nice_orthogroups_file = 'orthogroups_nice.txt'


with open(orthomcl_config) as f:
    conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')
    ortholog_table = conf['orthologTable'].strip()
    in_paralog_table = conf['inParalogTable'].strip()
    coortholog_table = conf['coOrthologTable'].strip()
    similar_sequeces_table = conf['similarSequencesTable'].strip()
    inter_taxon_match_view = conf['interTaxonMatchView'].strip()
    best_hit_table = 'BestHit'
    best_hit_taxon_score_table = 'BestQueryTaxonScore'


def find_genes(assembly):
    return Step(
        'Finding genes',
         run=cmdline('prodigal',
             parameters=[
                 '-a', predicted_proteins,
                 '-o', 'intermediate/predicted.gbk',
                 '-i', assembly]),
         req_files=assembly,
         prod_files=[predicted_proteins])

def step_fetching_annotations_for_species(specied_list, proxy):
    return Step(
        'Fetching annotations',
         run=lambda: fetch_annotations_for_species_from_ftp(
             annotations_dir, specied_list, proxy),
         prod_files=[annotations_dir])

def step_fetch_annotations_for_ids(ids_list):
    return Step(
        'Fetching annotations',
         run=lambda: fetch_annotations_for_ids(annotations_dir, ids_list),
         prod_files=[annotations_dir])

def step_make_proteomes(annotations=None):
    return Step(
        'Preparing proteomes',
         run=lambda: make_proteomes(annotations or annotations_dir, proteomes_dir),
         prod_files=[proteomes_dir])

def step_adjust_proteomes(proteomes_files, id_field=1):
    return Step(
        'Adjusting proteomes',
         run=lambda: adjust_proteomes(proteomes_files, proteomes_dir,
                                      id_field),
         req_files=proteomes_files,
         prod_files=[proteomes_dir])

def filter_proteomes():
    return Step(
        'Filtering proteomes',
         run=cmdline(join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
             parameters=[proteomes_dir,
                         10, 20,
                         good_proteins,
                         poor_proteins]),
         req_files=[proteomes_dir],
         prod_files=[good_proteins,
                     poor_proteins])

def make_blast_db():
    return Step(
        'Making blast database',
         run=cmdline(
             'makeblastdb',
              parameters=[
                '-in', good_proteins,
                '-input_type', 'fasta',
                '-out', blast_db,
                '-dbtype', 'prot'],
              stdout='log'),
         req_files=[good_proteins],
         prod_files=[blast_db + '.' + ext for ext in ['phr', 'pin', 'psq']])

def blast(threads):
    parameters = [
        '-query', good_proteins,
        '-db', blast_db,
        '-out', blast_out,
        '-outfmt', 6,  # tabular
        '-seg', 'yes',
        '-soft_masking', 'true',
        '-evalue', 1e-5,
        '-dbsize', BLAST_DBSIZE]

    def run():
        res = cmdline('blastp',
                      parameters + ['-num_threads', threads],
                      stdout=None)()
        if res == -6:
            log.info('')
            log.warn('Warning: blast refused to run multithreaded, running single-threaded instead.')
            res = cmdline('blastp', parameters)()
        return res

    return Step(
        'Blasting',
         run=run,
         req_files=[good_proteins],
         prod_files=[blast_out])

def parse_blast_restults():
    return Step(
        'Parsing blast results',
         run=cmdline(join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
            parameters=[blast_out, proteomes_dir],
            stdout=similar_sequences),
         req_files=[proteomes_dir, blast_out],
         prod_files=[similar_sequences])

def clean_database(suffix):
    return Step(
        'Cleaning database',
         run=lambda: clean_db(suffix))

def install_schema(suffix):
    return Step(
        'Installing schema',
         run=cmdline(join(orthomcl_bin_dir, 'orthomclInstallSchema.pl'),
                     parameters=[
                         orthomcl_config,
                         sql_log,
                         suffix],
                     stderr='log'),
         req_files=[orthomcl_config],
         prod_tables=[
             ortholog_table + suffix,
             in_paralog_table + suffix,
             coortholog_table + suffix,
             similar_sequeces_table + suffix,
             inter_taxon_match_view + suffix])

def load_blast_results(suffix):
    return Step(
        'Loading blast results into the database',
         run=cmdline(join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
                     parameters=[
                         orthomcl_config,
                         similar_sequences,
                         suffix],
                     stderr='log'),
         req_files=[orthomcl_config,
                    similar_sequences],  # and initialized database
         prod_files=[])  # loads blast results into the db)

def find_pairs(suffix):
    return Step(
        'Finding pairs',
         run=cmdline(
             join(orthomcl_bin_dir, 'orthomclPairs.pl'),
             parameters=[orthomcl_config,
                         pairs_log,
                         'cleanup=yes',
                         'suffix=' + suffix],
             stderr='log'),
         req_files=[orthomcl_config],
         req_tables=[in_paralog_table + suffix,
                     ortholog_table + suffix,
                     coortholog_table + suffix],
         prod_files=[])  # populates InParalog, Ortholog, CoOrtholog)

def dump_pairs_to_files(suffix):
    return Step(
        'Dump pairs files',
         run=cmdline(join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
             parameters=[orthomcl_config,
                         relpath(mcl_input, intermediate_dir),
                         intermediate_dir,
                         suffix],
             stderr='log'),
         req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
         req_tables=[in_paralog_table + suffix,
                     ortholog_table + suffix,
                     coortholog_table + suffix],
         prod_files=[mcl_input,
                     pairs_dir,
                     pairs_orthologs,
                     pairs_inparalogs,
                     pairs_coorthologs])

def mcl(inflation=1.5):
    return Step(
        'MCL',
         run=cmdline('mcl',
             parameters=[mcl_input,
                         '--abc',
                         '-I', str(inflation),
                         '-o', mcl_output],
             stderr='log',
             stdout='log'),
         req_files=[mcl_input],
         prod_files=[mcl_output])

def step_save_orthogroups(annotations=None):
    return Step(
        'Saving orthogroups',
        run=lambda: save_orthogroups(
            annotations or annotations_dir, mcl_output, orthogroups_file, nice_orthogroups_file),
        req_files=[mcl_output],
        prod_files=[orthogroups_file, nice_orthogroups_file])

#def groups_to_files(prefix, start_id):
#    return Step(
#        'MCL groups to files',
#         cmd=join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
#         req_files=[mcl_output],
#         prod_files=[groups_file],
#         parameters=[prefix + '_', start_id],
#         stdin=mcl_output,
#         stdout=groups_file)
#
#def signletones_to_files():
#    return Step(
#        'MCL singletones to files',
#         cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
#         req_files=[good_proteins, groups_file],
#         prod_files=[singletons_file],
#         parameters=[good_proteins, groups_file],
#         stdout=singletons_file)