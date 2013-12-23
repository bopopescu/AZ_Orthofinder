#!/usr/bin/env python

import sys
import logging
from random import randint
from shutil import copy
from os import chdir, mkdir, remove, access, X_OK, environ, getcwd
from os.path import join, isfile, exists, basename, isdir, relpath, split, pathsep, dirname, realpath
from clean_db import clean_db

import src.config as config
from src.fetch_annotations import fetch_annotations_species_name_ftp, fetch_annotations_ids
from src.make_proteomes import make_proteomes
from src.Workflow import Step, Workflow

log = logging.getLogger(config.log_fname)
orthomcl_config = config.orthomcl_config
orthomcl_bin_dir = config.orthomcl_bin_dir
filepath = dirname(realpath(__file__))


def make_workflow_id__from_species_names(species_names=None):
    if not species_names:
        return str(randint(1000, 9999))

    sn_words = ' '.join(species_names).split()
    if len(sn_words) == 1:
        return sn_words[0].replace('-', '')[:4]
    if len(sn_words) >= 2:
        return ''.join(w[0] for w in sn_words)


def make_workflow_id(working_dir=None):
    if not working_dir:
        return str(randint(1000, 9999))

    return basename(working_dir)[:4]


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous.')
    op.add_argument('--species-file', dest='species_file',
                    help='File with a list of organism names as in Genbank. '
                         'For example, "Salmonella enterica subsp. enterica '
                         'serovar Typhi str. P-stx-12".'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('--ids-file', dest='ids_file',
                    help='File with reference ids to fetch from Genbank.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('--annotations-dir', dest='annotations_dir',
                    help='Directory with .gb files.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('-o', dest='out_dir', required=True,
                    help='The directory that will contain the resulting group.txt file, '
                         'as well as intermediate results.')

    #op.add_argument('-w', '--overwrite', dest='overwrite', action='store_true', default=False,
    #                help='By default, the tool reuses existing intermediate results.'
    #                     'This option makes the tool overwrite any existing data.')

    op.add_argument('--start-from', dest='start_from', default=0,
                    help='Start from the specified step. Either name (see log.txt) or "uselog".'
                         'If "uselog", the last "Done" record in log.txt will be searched.')

    op.add_argument('-t', dest='threads', default=1,
                    help='Number of threads to run Blast.')

    op.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)

    op.add_argument('--ask', '--ask-each-step',
                    dest='ask_each_step', action='store_true', default=False,
                    help='Wait for user to press ke every time before proceed to next step.')

    op.add_argument('--proxy', dest='proxy', default=None, help='Proxy for FTP, for example: '
                                                                '--proxy 198.260.1.1:3333')
    usage = 'find_orthologs [--species-file FILE] \n' \
            '                      [--ids-file FILE] \n' \
            '                      [--annotations-dir DIR] \n' \
            '                       -o OUTPUT_DIRECTORY \n' \
            '                      [--start-from STEP_NAME]'
    op.usage = usage

    params = op.parse_args(args)

    if not params.species_file and not params.ids_file and not params.annotations_dir \
            and not params.start_from:
        print >> sys.stderr, 'Either --species-file, --ids-file, --annotations-dir, ' + \
                             'or --start-from has to be specified.'
        op.exit(1)

    if params.species_file and not isfile(params.species_file):
        print >> sys.stderr, 'Species list file ' + params.species_file + ' does not exist or is directory.'
        exit(1)

    if params.ids_file and not isfile(params.ids_file):
        print >> sys.stderr, 'Reference ids file ' + params.ids_file + ' does not exist or is directory.'
        exit(1)

    if params.annotations_dir and not isdir(params.annotations_dir):
        print >> sys.stderr, 'Annotations directory ' + params.annotations_dir + ' does not exist or is file.'
        exit(1)

    if params.start_from == 'uselog':
        if not isfile(join(params.out_dir, config.log_fname)):
            print >> sys.stderr, 'No %s in %s. Either check your path, or ' \
                                 'change the --start-from option' % \
                                 (config.log_fname, params.out_dir)
            exit(1)

    return params


def which(program):
    def is_exe(fpath):
        return isfile(fpath) and access(fpath, X_OK)

    fpath, fname = split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in environ["PATH"].split(pathsep):
            path = path.strip('"')
            exe_file = join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def run_workflow(working_dir, overwrite,
                 specied_list, ids_list, user_annotations_dir,
                 ask_before=False,
                 start_after=None, start_from=None, threads=1,
                 proxy=None,
                 inflation=1.5, first_id=1000):

    if not which('blastp') or not which('mcl'):
        if not which('blastp'):
            log.error('blastp installation required.')
        if not which('mcl'):
            log.error('mcl installation required.')
        return 3

    workflow_id = make_workflow_id(working_dir)
    log.info('Workflow id is "' + workflow_id + '"')
    log.info('')

    if not exists('intermediate'):
        mkdir('intermediate')

    #clean_db_script = join(filepath, 'clean_db.py')

    proteomes_dir       = 'proteomes'
    annotations_dir     = 'annotations'
    sql_log             = 'intermediate/log.sql'
    good_proteins       = 'intermediate/good_proteins.fasta'
    poor_proteins       = 'intermediate/poor_proteins.fasta'
    blast_db            = 'intermediate/blastdb'
    blast_out           = 'intermediate/blasted.tsv'
    similar_sequences   = 'intermediate/similar_sequences.txt'
    pairs_log           = 'intermediate/orthomclpairs.log'
    mcl_input           = 'intermediate/mcl_input'
    mcl_output          = 'intermediate/mcl_output'
    groups_file         = 'groups.txt'
    singletons_file     = 'singletons.txt'

    with open(orthomcl_config) as f:
        conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')
        ortholog_table = conf['orthologTable'].strip() + '_' + workflow_id
        in_paralog_table = conf['inParalogTable'].strip() + '_' + workflow_id
        coortholog_table = conf['coOrthologTable'].strip() + '_' + workflow_id
        similar_sequeces_table = conf['similarSequencesTable'].strip() + '_' + workflow_id
        inter_taxon_match_view = conf['interTaxonMatchView'].strip() + '_' + workflow_id
        best_hit_table = 'BestHit' + '_' + workflow_id
        best_hit_taxon_score_table = 'BestQueryTaxonScore' + '_' + workflow_id

    steps = []

    if specied_list:
        log.debug('Using species list: ' + str(specied_list))
        steps.append(
            Step('Fetching annotations',
                 cmd=fetch_annotations_species_name_ftp,
                 prod_files=[annotations_dir],
                 parameters=[annotations_dir, specied_list, proxy]))

    elif ids_list:
        log.debug('Using ref ids: ' + str(ids_list))
        steps.append(
            Step('Fetching annotations',
                 cmd=fetch_annotations_ids,
                 prod_files=[annotations_dir],
                 parameters=[annotations_dir, ids_list]))

    elif user_annotations_dir:
        log.debug('Using user_annotations_dir: ' + user_annotations_dir)
        annotations_dir = user_annotations_dir

    elif start_from == 0:
        log.error('Either species names, reference ids, annotations directory, '
                  'or step to start from has to be be speciefied.')
        exit(1)

    steps.extend([
        Step('Preparing proteins',
             cmd=make_proteomes,
             req_files=[annotations_dir],
             prod_files=[proteomes_dir],
             parameters=[
                annotations_dir,
                workflow_id,
                proteomes_dir]),

        Step('Filtering fasta',
             cmd=join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
             req_files=[proteomes_dir],
             prod_files=[good_proteins, poor_proteins],
             parameters=[
                proteomes_dir,
                10,
                20,
                good_proteins,
                poor_proteins]),

        Step('Making blast database',
             cmd='makeblastdb',
             req_files=[good_proteins],
             prod_files=[blast_db + '.' + ext for ext in ['phr', 'pin', 'psq']],
             parameters=[
                '-in', good_proteins,
                '-input_type', 'fasta',
                '-out', blast_db,
                '-dbtype', 'prot',
                '-title', basename(working_dir)],
             stdout='log'),

        Step('Blasting all vs. all',
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
                '-dbsize', config.BLAST_DBSIZE]),

        Step('Parsing blast results',
             cmd=join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
             req_files=[proteomes_dir, blast_out],
             prod_files=[similar_sequences],
             parameters=[
                blast_out,
                proteomes_dir],
             stdout=similar_sequences),

        Step('Cleaning database',
             cmd=clean_db,
             parameters=[workflow_id])
        if overwrite else None,

        Step('Installing schema',
             cmd=join(orthomcl_bin_dir, 'orthomclInstallSchema.pl'),
             req_files=[orthomcl_config],
             prod_files=[],  # drops and creates SimilarSequences, InParalog,
                             # Ortholog, CoOrtholog, InterTaxonMatch view
             prod_tables=[
                ortholog_table,
                in_paralog_table,
                coortholog_table,
                similar_sequeces_table,
                inter_taxon_match_view],
             parameters=[
                orthomcl_config,
                sql_log,
                ('_' + workflow_id if workflow_id else '')],
             stderr='log'),

        Step('Loading blast results into the database',
             cmd=join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
             req_files=[orthomcl_config, similar_sequences],  # and initialized database
             prod_files=[],  # loads blast results into the db
             parameters=[
                orthomcl_config,
                similar_sequences,
                ('_' + workflow_id) if workflow_id else ''],
             stderr='log'),

        Step('Finding pairs',
             cmd=join(orthomcl_bin_dir, 'orthomclPairs.pl'),
             req_files=[orthomcl_config],
             req_tables=[in_paralog_table, ortholog_table, coortholog_table],
             prod_files=[],  # populates InParalog, Ortholog, CoOrtholog
             parameters=[
                orthomcl_config,
                pairs_log,
                'cleanup=no',
                #'startAfter=useLog' if not overwrite else '',
                ('suffix=_' + workflow_id) if workflow_id else ''],
             stderr='log'),

        Step('Dump pairs files',
             cmd=join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
             req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
             req_tables=[in_paralog_table, ortholog_table, coortholog_table],
             prod_files=[
                mcl_input,
                'pairs',
                'pairs/potentialOrthologs.txt',
                'pairs/potentialInparalogs.txt',
                'pairs/potentialCoorthologs.txt'],
             parameters=[
                orthomcl_config,
                mcl_input,
                working_dir,  # current dir. no effect, just to use next positional argument
                ('_' + workflow_id) if workflow_id else ''],
             stderr='log'),

        Step('MCL',
             cmd='mcl',
             req_files=[mcl_input],
             prod_files=[mcl_output],
             parameters=[
                mcl_input,
                '--abc',
                '-I', str(inflation),
                '-o', mcl_output],
             stderr='log',
             stdout='log'),

        Step('MCL groups to files',
             cmd=join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
             req_files=[mcl_output],
             prod_files=[groups_file],
             parameters=[
                workflow_id + '_',
                str(first_id)],
             stdin=mcl_output,
             stdout=groups_file),

        Step('MCL singletones to files',
             cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
             req_files=[good_proteins, groups_file],
             prod_files=[singletons_file],
             parameters=[
                good_proteins,
                groups_file],
             stdout=singletons_file),
    ])

    workflow = Workflow(steps)
    result = workflow.run(start_after, start_from, overwrite, ask_before)
    if result == 0:
        log.info('Done.')
        log.info('Log in ' + join(working_dir, config.log_fname))
        log.info('Groups in ' + join(working_dir, groups_file))
    return result


def set_up_config():
    with open(config.config) as cf:
        conf = dict(l.strip().split('=', 1) for l
                    in cf.readlines() if l.strip()[0] != '#')
        log.debug('Read conf: ' + str(conf))
    with open(orthomcl_config) as ocf:
        omcl_conf = dict(l.strip().split('=', 1) for l
                         in ocf.readlines() if l.strip()[0] != '#')
    omcl_conf['dbConnectString'] = \
        'dbi:mysql:database=orthomcl;host=127.0.0.1;port=%s;mysql_local_infile=1' % \
            conf['db_port']
    omcl_conf['dbLogin'] = conf['db_login']
    omcl_conf['dbPassword'] = conf['db_password']

    with open(orthomcl_config, 'w') as ocf:
        ocf.writelines('='.join(item) + '\n' for item in omcl_conf.items())


def read_list(file, out_dir):
    if not file:
        return None
    with open(file) as f:
        results = [l.strip() for l in f.readlines() if l.strip()]
    if isfile(join(out_dir, basename(file))):
        remove(join(out_dir, basename(file)))
    copy(file, out_dir)
    return results


def main(args):
    params = parse_args(args)

    if not isdir(params.out_dir):
        mkdir(params.out_dir)

    config.set_up_logging(params.debug, params.out_dir)

    species_list = read_list(params.species_file, params.out_dir)
    ref_id_list = read_list(params.ids_file, params.out_dir)
    if params.annotations_dir:
        params.annotations_dir = join(getcwd(), params.annotations_dir)

    start_after = None
    start_from = params.start_from
    if params.start_from == 'uselog':
        with open(join(params.out_dir, config.log_fname)) as f:
            for l in f.readlines():
                if 'Done' in l:
                    start_after = l[41:].strip()
                    start_from = None

    set_up_config()

    log.info('Changing to %s' % params.out_dir)
    chdir(params.out_dir)

    run_workflow(working_dir=params.out_dir,
                 #overwrite=params.overwrite or start_after != 0 or start_from != 0,
                 overwrite=True,
                 specied_list=species_list,
                 ids_list=ref_id_list,
                 user_annotations_dir=params.annotations_dir,
                 ask_before=params.ask_each_step,
                 start_after=start_after,
                 start_from=start_from,
                 threads=params.threads,
                 proxy=params.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])