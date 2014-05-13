from genericpath import isfile
from os import remove, listdir
from os.path import basename, join, relpath, exists, isdir, realpath
from shutil import rmtree

from Workflow import Step, cmdline
from src.db_connection import DbCursor
from utils import check_install_mcl, check_installed_tools
from process_assembly import filter_assembly
from save_orthogroups import save_orthogroups
from make_proteomes import make_proteomes, adjust_proteomes
from fetch_annotations import fetch_annotations_for_species_from_ftp, fetch_annotations_for_ids
from config import orthomcl_config_final_path, orthomcl_bin_dir, BLAST_DBSIZE
import config
from clean_db import clean_db

import config
import logging
log = logging.getLogger(config.log_fname)


def check_results_existence():
    existing = [obj for obj in [
        config.proteomes_dir,
        config.annotations_dir,
        config.intermediate_dir,
        config.groups_file,
        config.orthogroups_file,
        config.nice_orthogroups_file] if exists(obj)]

    if existing:
        log.warn('The directory contains previous results. Do you want to overwrite it? '
                 '(You can run with the --overwrite option to avoid this warning.)')
        raw_input('Press any key to overwrite and continue, or Ctrl-C to interrupt.\n> ')

    for obj in existing:
        if isfile(obj):
            remove(obj)
        if isdir(obj):
            rmtree(obj)


ortholog_table = 'Ortholog'
in_paralog_table = 'InParalog'
coortholog_table = 'CoOrtholog'
similar_sequeces_table = 'SimilarSequences'
inter_taxon_match_view = 'interTaxonMatchView'
best_hit_table = 'BestHit'
best_hit_taxon_score_table = 'BestQueryTaxonScore'


#def step_fetching_annotations_for_species(specied_list, proxy):
#    return Step(
#        'Fetching annotations',
#         run=lambda: fetch_annotations_for_species_from_ftp(
#             annotations_dir, specied_list, proxy),
#         prod_files=[annotations_dir])
#
#def step_fetch_annotations_for_ids(ids_list):
#    return Step(
#        'Fetching annotations',
#         run=lambda: fetch_annotations_for_ids(annotations_dir, ids_list),
#         prod_files=[annotations_dir])
#
#def step_make_proteomes(annotations=None):
#    return Step(
#        'Preparing proteomes',
#         run=lambda: make_proteomes(annotations or annotations_dir, proteomes_dir),
#         prod_files=[proteomes_dir])


#def filter_new_proteome(assembly_name, min_length=10, max_percent_stop=20):
#    return Step(
#        'Filtering proteome',
#         run=cmdline(join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
#             parameters=[join(proteomes_dir, assembly_name + '.faa'),
#                         min_length, max_percent_stop,
#                         good_proteins,
#                         poor_proteins]),
#         req_files=[proteomes_dir, join(proteomes_dir, assembly_name + '.faa')],
#         prod_files=[good_proteins, poor_proteins])


######################################################################
#def step_adjust_proteomes(proteomes_files, id_field=1):
#    return Step(
#        'Adjusting proteomes',
#         run=lambda: adjust_proteomes(proteomes_files, proteomes_dir,
#                                      id_field),
#         req_files=proteomes_files,
#         prod_files=[proteomes_dir])

def filter_proteomes(min_length=10, max_percent_stop=20):
    return Step(
        'Filtering proteomes: min length = ' + str(min_length) +
        ', max percent of stop codons = ' + str(max_percent_stop),
        run=cmdline(
            'perl' + join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
            parameters=[
                realpath(config.proteomes_dir),
                min_length, max_percent_stop,
                realpath(config.good_proteins),
                realpath(config.poor_proteins)]),
        req_files=[config.proteomes_dir],
        prod_files=[config.good_proteins, config.poor_proteins])

def filter_proteomes_split(min_length=10, max_percent_stop=20, jobs=30):
    def proc():
        # 1. iterate each fasta, count total number of good proteins
        # 2. devide by N number of jobs, get K
        # 3. iterate each fasta, save N good_protein files with K proteins
        # 4. return multiple good_proteins
        pass

    return Step(
        'Filtering proteomes: min length = ' + str(min_length) +
        ', max percent of stop codons = ' + str(max_percent_stop),
        run=proc,
        req_files=[config.proteomes_dir],
        prod_files=[config.good_proteins, config.poor_proteins])

def make_blast_db():
    return Step(
        'Making blast database',
        run=cmdline(
            'makeblastdb',
             parameters=[
               '-in', realpath(config.good_proteins),
               '-input_type', 'fasta',
               '-out', realpath(config.blast_db),
               '-dbtype', 'prot'],
             stdout='log',
             stderr='log'),
        req_files=[config.good_proteins],
        prod_files=[config.blast_db + '.' + ext for ext in ['phr', 'pin', 'psq']])

def blast(threads, new_good_proteomes=None, evalue=1e-5):
    if new_good_proteomes:
        parameters = [
            '-query', realpath(new_good_proteomes),
            '-db', realpath(config.blast_db),
            '-out', realpath(config.blast_out + '_2'),
            '-outfmt', 6,  # tabular
            '-seg', 'yes',
            '-soft_masking', 'true',
            '-evalue', evalue,
            '-dbsize', BLAST_DBSIZE]
    else:
        parameters = [
            '-query', realpath(config.good_proteins),
            '-db', realpath(config.blast_db),
            '-out', realpath(config.blast_out),
            '-outfmt', 6,  # tabular
            '-seg', 'yes',
            '-soft_masking', 'true',
            '-evalue', evalue,
            '-dbsize', BLAST_DBSIZE]

    def run():
        if not check_installed_tools(['blastp'], only_warn=False):
            return 1

        res = cmdline(
            'blastp', parameters + ['-num_threads', threads],
            ignore_lines_by_pattern=r'.* at position .* replaced by .*')()
        if res == -6:
            log.warn('')
            log.warn('   WARNING: blast refused to run multithreaded, '
                     'running single-threaded instead.')
            res = cmdline('blastp', parameters)()

        if new_good_proteomes:
            log.info('   Appending ' + config.blast_out + '_2 to ' + config.blast_out)

            with open(config.blast_out, 'a') as b_out:
                with open(config.blast_out + '_2') as b_out_2:
                    b_out.write(b_out_2.read())
        return res

    return Step(
        'Blasting',
        run=run,
        req_files=[config.good_proteins])

def parse_blast_results():
    return Step(
        'Parsing blast results',
        run=cmdline(
            join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
            parameters=[realpath(config.blast_out), realpath(config.proteomes_dir)],
            stdout=realpath(config.similar_sequences)),
        req_files=[config.proteomes_dir, config.blast_out],
        prod_files=[config.similar_sequences])

def clean_database(suffix):
    return Step(
        'Cleaning database',
        run=lambda: clean_db(suffix))

def install_schema(suffix):
    return Step(
        'Installing schema',
        run=cmdline(
            join(orthomcl_bin_dir, 'orthomclInstallSchema.pl'),
            parameters=[
                realpath(orthomcl_config_final_path),
                realpath(config.sql_log),
                suffix],
            stderr='log'),
        req_files=[orthomcl_config_final_path],
        prod_tables=[
            ortholog_table + suffix,
            in_paralog_table + suffix,
            coortholog_table + suffix,
            similar_sequeces_table + suffix,
            inter_taxon_match_view + suffix])

def load_blast_results(suffix):
    def run():
        with DbCursor() as cursor:
            for tbl in [
                similar_sequeces_table + suffix,
            ]:
                try:
                    log.debug('   Cleaning the %s table.' % tbl)
                    try:
                        cursor.execute('select 1 from %s limit 1;' % tbl)
                    except:
                        pass
                    log.debug('   select 1 from ' + tbl + ' limit 1; '
                              'result=' + str(cursor.fetchone()))
                    try:
                        cursor.execute('delete from %s;' % tbl)
                        cursor.execute('select 1 from %s limit 1;' % tbl)
                    except:
                        pass
                    log.debug('   select 1 from ' + tbl + ' limit 1; '
                              'result=' + str(cursor.fetchone()))
                    log.debug('')

                except Exception, e:
                    log.exception(e)

        return cmdline(join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
            parameters=[
                realpath(orthomcl_config_final_path),
                realpath(config.similar_sequences),
                suffix,
                ])()

    return Step(
        'Loading blast results into the database',
        run=run,
        req_files=[orthomcl_config_final_path,
                   config.similar_sequences],  # and initialized database
        prod_files=[])  # loads blast results into the db)

def find_pairs(suffix):
    def run():
        with DbCursor() as cursor:
            for tbl in [
                in_paralog_table + suffix,
                ortholog_table + suffix,
                coortholog_table + suffix,
            ]:
                try:
                    log.debug('   Cleaning the %s table.' % tbl)
                    try:
                        cursor.execute('select 1 from %s limit 1;' % tbl)
                    except:
                        pass
                    log.debug('   select 1 from ' + tbl + ' limit 1; '
                              'result=' + str(cursor.fetchone()))
                    try:
                        cursor.execute('delete from %s;' % tbl)
                        cursor.execute('select 1 from %s limit 1;' % tbl)
                    except:
                        pass
                    log.debug('   select 1 from ' + tbl + ' limit 1; '
                              'result=' + str(cursor.fetchone()))
                    log.debug('')

                except Exception, e:
                    log.exception(e)

        #res = cmdline(
        #    join(orthomcl_bin_dir, 'orthomclPairs.pl'),
        #    parameters=[
        #        orthomcl_config,
        #        config.pairs_log,
        #        'cleanup=only',
        #        'suffix=' + (suffix if suffix else '*')])()
        #log.info('   Cleaning: ' + str(res))

        return cmdline(
            join(orthomcl_bin_dir, 'orthomclPairs.pl'),
            parameters=[
                realpath(orthomcl_config_final_path),
                realpath(config.pairs_log),
                'cleanup=no',
                'suffix=' + (suffix if suffix else '*')])()

    return Step(
        'Finding pairs',
        run=run,
        req_files=[orthomcl_config_final_path],
        req_tables=[
            in_paralog_table + suffix,
            ortholog_table + suffix,
            coortholog_table + suffix,
            similar_sequeces_table + suffix],
        prod_tables=[tbl_name + suffix for tbl_name in [
            'BestHit',
            'BestQueryTaxonScore',
            'BestInterTaxonScore',
            'BetterHit',
            'CoOrthNotOrtholog',
            'CoOrthologTaxon',
            'CoOrthologCandidate',
            'CoOrthologAvgScore',
            'CoOrthologTemp',
            'InParalog2Way',
            'InParalogAvgScore',
            'InParalogTemp',
            'InParalogTaxonAvg',
            'InParalogOrtholog',
            'InplgOrthTaxonAvg',
            'InplgOrthoInplg',
            'OrthologAvgScore',
            'OrthologTemp',
            'Ortholog2Way',
            'OrthologTaxon',
            'OrthologUniqueId',
            'UniqSimSeqsQueryId',
        ]],
        prod_files=[])  # populates InParalog, Ortholog, CoOrtholog)

def dump_pairs_to_files(suffix):
    def run():
        res = cmdline(
            join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
             parameters=[realpath(orthomcl_config_final_path),
                         realpath(config.mcl_input),
                         realpath(config.intermediate_dir),
                         suffix],
             stderr='log')()

        with DbCursor() as cursor:
            for tbl in [
                in_paralog_table + suffix,
                ortholog_table + suffix,
                coortholog_table + suffix,
            ]:
                try:
                    log.debug('   Cleaning the %s table.' % tbl)
                    cursor.execute('select 1 from %s limit 1;' % tbl)
                    log.debug('   ' + str(cursor.fetchone()))
                    cursor.execute('delete from %s;' % tbl)
                    cursor.execute('select 1 from %s limit 1;' % tbl)
                    log.debug('   ' + str(cursor.fetchone()))
                    log.debug('')

                except Exception, e:
                    log.exception(e)
        return res

    return Step(
        'Dumping pairs files',
        run=run,
        req_files=[orthomcl_config_final_path],  # and populated InParalog, Ortholog, CoOrtholog tables
        #req_tables=[in_paralog_table + suffix,
        #            ortholog_table + suffix,
        #            coortholog_table + suffix],
        prod_files=[config.mcl_input,
                    config.pairs_dir,
                    config.pairs_orthologs,
                    config.pairs_inparalogs,
                    config.pairs_coorthologs])

def mcl(debug, inflation=1.5):
    def run():
        mcl_bin_path, res = check_install_mcl(debug, only_warn=False)
        if mcl_bin_path is None:
            return res

        return cmdline(
            mcl_bin_path,
            parameters=[
                realpath(config.mcl_input),
                '--abc',
                '-I', str(inflation),
                '-o', realpath(config.mcl_output)],
            start_ignoring_from=r'Please cite:.*',
            stderr='log',
            stdout='log')()

    return Step(
        'MCL',
        run=run,
        req_files=[config.mcl_input],
        prod_files=[config.mcl_output])

def step_save_orthogroups(added_proteomes_dir=None,
                          annotations=None, internet_on=True):
    def run():
        if added_proteomes_dir:
            added_proteomes_files = [
                join(added_proteomes_dir, prot)
                for prot in listdir(added_proteomes_dir)
                if prot and prot[0] != '.']
        else:
            added_proteomes_files = []

        return save_orthogroups(
            added_proteomes_files, annotations or config.annotations_dir, config.mcl_output,
            config.orthogroups_file, config.nice_orthogroups_file,
            config.short_orthogroups_file, config.assembly_singletones_file,
            config.singletone_dir)

    prod_files = [
        config.orthogroups_file,
        config.nice_orthogroups_file,
        config.short_orthogroups_file,
        config.assembly_singletones_file]

    return Step(
       'Saving orthogroups',
       run=run,
       req_files=[config.mcl_output],
       prod_files=prod_files)

def groups_to_files(prefix, start_id=0):
    return Step(
        'MCL groups to files',
        run=cmdline(
            join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
            parameters=[prefix + '_', start_id],
            stdin=config.mcl_output,
            stdout=config.groups_file),
        req_files=[config.mcl_output],
        prod_files=[config.groups_file])

#def signletones_to_files():
#    return Step(
#        'MCL singletones to files',
#         cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
#         req_files=[good_proteins, groups_file],
#         prod_files=[singletons_file],
#         parameters=[good_proteins, groups_file],
#         stdout=singletons_file)