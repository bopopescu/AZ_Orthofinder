from random import randint
from shutil import rmtree, copyfile, copy
from subprocess import call
from os import chdir, getcwd, listdir, mkdir, remove
from os.path import join, realpath, isfile, exists, dirname, basename, normpath, isdir
from itertools import ifilterfalse, ifilter, izip, count
from mysql.connector import errorcode
import mysql.connector
import sys

from src.fetch_annotations import fetch_annotations
from src.make_proteomes import make_proteomes
from src.db_connection import DbCursor

import logging
log = logging.getLogger('orthofinder')


orthomcl_config = join(dirname(realpath(__file__)), 'src/orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'src/orthomcl_software/bin')

log_file = 'log.txt'


def set_up_logging(debug, working_dir):
    logger = logging.getLogger('orthofinder')
    logger.setLevel(logging.DEBUG)

    class InfoFilter(logging.Filter):
        def filter(self, rec):
            return rec.levelno in (logging.DEBUG, logging.INFO)

    console_formatter = logging.Formatter(
        '%(asctime)-15s  %(message)s' if debug else '%(message)s',
        datefmt='%c')

    std = logging.StreamHandler(sys.stdout)
    std.setLevel(logging.DEBUG if debug else logging.INFO)
    std.addFilter(InfoFilter())
    std.setFormatter(console_formatter)
    logger.addHandler(std)

    err = logging.StreamHandler(sys.stderr)
    err.setLevel(logging.WARN)
    err.setFormatter(console_formatter)
    logger.addHandler(err)

    fh = logging.FileHandler(join(working_dir, log_file), 'a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)-15s  %(levelname)-8s  %(message)s',
        datefmt='%c'))
    logger.addHandler(fh)


def make_workflow_id(species_names=None):
    if not species_names:
        return randint(1000, 9999)

    sn_words = ' '.join(species_names).split()
    if len(sn_words) == 1:
        return sn_words[0].replace('-', '')[:11]
    if len(sn_words) >= 2:
        return ''.join(w[0] for w in sn_words)


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous.')
    op.add_argument('--species-file', dest='species_file',
                    help='File with a list of full organism names as in Genbank. '
                         'For example, "Salmonella enterica subsp. enterica '
                         'serovar Typhi str. P-stx-12".')

    op.add_argument('--ids', dest='ids', help='Comma-separated list of reference ids.')

    op.add_argument('-o', dest='out_dir',
                    help='The directory that will contain the resulting group.txt file,'
                         'and the intermediate results in data.')

    op.add_argument('--overwrite', dest='overwrite', action='store_true', default=False,
                    help='By default, the tool reuses existing intermediate results.'
                         'This option makes the tool preliminarily remove existing data')

    op.add_argument('--start-from', dest='start_from', default=0,
                    help='Start from this step. Either name or "uselog".'
                         'If "uselog", find the last "Done" record in log.txt file.')

    op.add_argument('-t', dest='threads', default=1,
                    help='Number of threads to run Blast.')

    op.add_argument('--blast-dbsize', dest='blast_dbsize', default=100000,
                    help='Maximum intended number of proteins to blast '
                         'after all incremental reruns.')

    op.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)

    op.add_argument('--ask-each-step',
                    dest='ask_each_step', action='store_true', default=False,
                    help='Ask user every time before proceed to next step.')

    op.add_argument('--proxy', dest='proxy', default=None, help='Set up proxy for FTP.')

    params = op.parse_args(args)
    if params.start_from == 'uselog':
        if not isfile(join(params.out_dir, log_file)):
            print >> sys.stderr, 'No %s in %s. Either check your path, or ' \
                                 'change the --start-from=uselog option' % (log_file, params.out_dir)
            exit(1)
    return params


class Step:
    def __init__(self, name, cmd,
                 req_files=None, prod_files=None,
                 req_tables=None, prod_tables=None,
                 parameters=None, stdin=None, stdout=None):
        self.name = name
        self.command = cmd
        self.req_files = req_files or []
        self.req_tables = req_tables or []
        self.prod_files = prod_files or []
        self.prod_tables = prod_tables or []
        self.parameters = parameters or []
        self.stdin = stdin
        self.stdout = stdout

    def program(self):
        return basename(self.command)

    def run(self, overwrite=False, step_by_step=False):
        assert self.command

        # Checking existence of produced tables and files
        missing_prod_files = list(ifilterfalse(exists, self.prod_files))

        missing_prod_tables = []
        existing_prod_tables = []
        if self.prod_tables:
            with DbCursor() as cursor:
                for table in self.prod_tables:
                    try:
                        query = 'select count(*) from %s;' % table
                        cursor.execute(query)
                    except mysql.connector.Error, err:
                        log.debug('err.errno == errorcode.ER_TABLE_EXISTS_ERROR: ' +
                                  str(err.errno == errorcode.ER_TABLE_EXISTS_ERROR))
                        log.debug(err.msg)
                        missing_prod_tables.append(table)
                    else:
                        existing_prod_tables.append(table)

        if not overwrite:
            if self.prod_files and not missing_prod_files:
                log.info('   All files to be produced already exist: ' +
                         ', '.join(self.prod_files))
            if self.prod_tables and not missing_prod_tables:
                log.info('   All tables to be installed already exist: ' +
                         ', '.join(self.prod_tables))
            if not missing_prod_files and not missing_prod_tables:
                log.info('   Skipping')
                return 0

        # Checking requirements
        missing_req_files = list(ifilterfalse(exists, self.req_files))
        if missing_req_files:
            log.error('   ' + self.name + ' requires files ' +
                      ', '.join(missing_req_files))
            return 1

        missing_req_tables = []
        if self.req_tables:
            with DbCursor() as cursor:
                for table in self.req_tables:
                    try:
                        cursor.execute('select count(*) from %s;' % table)
                    except mysql.connector.Error, err:
                        missing_req_tables.append(table)
        if missing_req_tables:
            log.error('   ' + self.name + ' requires tables ' +
                      ', '.join(missing_req_tables) + ' installed')
            return 1

        # Removing existing data of overwrite
        existing_prod_files = list(ifilter(exists, self.prod_files))
        if overwrite and existing_prod_files:
            log.info('   Overwriting existing ' + ', '.join(existing_prod_files))
            for file in existing_prod_files:
                if isfile(file):
                    remove(file)
                if isdir(file):
                    rmtree(file)

        if overwrite and existing_prod_tables:
            with DbCursor() as cursor:
                for table in existing_prod_files:
                    try:
                        cursor.execute('drop table %s;' % table)
                    except mysql.connector.Error, err:
                        log.critical(err)

        # Running
        if step_by_step:
            raw_input('   Proceed?')

        if hasattr(self.command, '__call__'):
            return self.command(*self.parameters)
        else:
            commandline = ' '.join([self.command] + map(str, self.parameters))
            stdin_f, stdout_f = None, None
            if self.stdin:
                commandline += ' < ' + self.stdin
                stdin_f = open(self.stdin)
            if self.stdout:
                stdout_f = open(self.stdout, 'w')
                commandline += ' > ' + self.stdout
            log.info('   ' + commandline)
            try:
                return call([self.command] + map(str, self.parameters),
                            stdin=stdin_f, stdout=stdout_f)
            except KeyboardInterrupt:
                return 1


def run_workflow(working_dir,
                 species_names,
                 overwrite,
                 blast_dbsize,
                 ask_before=False,
                 start_after=None,
                 start_from=None,
                 threads=1,
                 proxy=None,
                 inflation=1.5,
                 first_id=1000):

    workflow_id = make_workflow_id(species_names)
    log.info('Workflow id is ' + workflow_id)
    log.info('')

    if not exists('intermediate'):
        mkdir('intermediate')

    proteomes_dir       = 'proteomes'
    annotations_dir     = 'annotations'
    sql_log             = 'intermediate/log.sql'
    good_proteins       = 'intermediate/good_proteins.fasta'
    poor_proteins       = 'intermediate/poor_proteins.fasta'
    blast_out           = 'intermediate/blasted.tsv'
    similar_sequences   = 'intermediate/similar_sequences.txt'
    pairs_log           = 'intermediate/orthomclpairs.log'
    mcl_input           = 'mclInput'
    mcl_output          = 'mclOutput'
    gene_pairs_dir      = 'gene_pairs'
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

    steps = [
        Step('Fetching annotations',
             cmd=fetch_annotations,
             prod_files=[annotations_dir],
             parameters=[
                annotations_dir,
                species_names,
                proxy,
                ]),

        Step('Making proteins',
             cmd=make_proteomes,
             req_files=[annotations_dir],
             prod_files=[proteomes_dir],
             parameters=[
                annotations_dir,
                species_names,
                proteomes_dir,
                ]),

        Step('Filtering fasta',
             cmd=join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
             req_files=[proteomes_dir],
             prod_files=[good_proteins, poor_proteins],
             parameters=[
                proteomes_dir,
                10,
                20,
                good_proteins,
                poor_proteins,
                ]),

        Step('Blasting all vs. all',
             cmd='blastp',
             req_files=[proteomes_dir, good_proteins],
             prod_files=[blast_out],
             parameters=[
                '-query', good_proteins,
                '-subject', good_proteins,
                '-out', blast_out,
                '-outfmt', 6,  # tabular
                '-evalue', 1e-5,
                '-num_descriptions', 10000,  # don't care value
                '-num_alignments', 10000,  # don't care value,
                '-num_threads', threads,
                '-dbsize', blast_dbsize,
                ]),

        Step('Parsing blast results',
             cmd=join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
             req_files=[proteomes_dir, blast_out],
             prod_files=[similar_sequences],
             parameters=[
                blast_out,
                proteomes_dir],
             stdout=similar_sequences),

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
                ('_' + workflow_id if workflow_id else ''),
                ]),

        Step('Loading blast results into the database',
             cmd=join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
             req_files=[orthomcl_config, similar_sequences],  # and initialized database
             prod_files=[],  # loads blast results into the db
             parameters=[
                orthomcl_config,
                similar_sequences,
                ('_' + workflow_id) if workflow_id else ''
                ]),

        Step('Finding pairs',
             cmd=join(orthomcl_bin_dir, 'orthomclPairs.pl'),
             req_files=[orthomcl_config],  # and initialized database
             prod_files=[],  # populates InParalog, Ortholog, CoOrtholog
             parameters=[
                orthomcl_config,
                pairs_log,
                'cleanup=no',
                #'startAfter=useLog' if not overwrite else '',
                ('suffix=_' + workflow_id) if workflow_id else '',
                ]),

        Step('Dump pairs files',
             cmd=join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
             req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
             prod_files=[
                mcl_input,
                'pairs',
                'pairs/potentialOrthologs.txt',
                'pairs/potentialInparalogs.txt',
                'pairs/potentialCoorthologs.txt'],
             parameters=[
                orthomcl_config,
                working_dir,  # current dir. no effect, just to use next positional argument
                ('_' + workflow_id) if workflow_id else '']),

        Step('MCL',
             cmd='mcl',
             req_files=[mcl_input],
             prod_files=[mcl_output],
             parameters=[
                mcl_input,
                '--abc',
                '-I', str(inflation),
                '-o', mcl_output,
                ]),

        Step('MCL groups to files',
             cmd=join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
             req_files=[mcl_output],
             prod_files=[groups_file],
             parameters=[
                workflow_id + '_',
                str(first_id)],
             stdin='mclOutput',
             stdout=groups_file),

        Step('MCL singletones to files',
             cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
             req_files=[good_proteins, groups_file],
             prod_files=[singletons_file],
             parameters=[
                good_proteins,
                groups_file],
             stdout=singletons_file),
    ]
    for i, step in izip(count(1), steps):
        if start_after is not None \
                and isinstance(start_after, basestring) \
                and start_after.lower() == step.name.lower():
            start_after = i
            continue

        if start_from is not None \
            and isinstance(start_from, basestring) \
            and start_from.lower() == step.name.lower():
            start_from = i

        if start_after is not None:
            if not isinstance(start_after, int) or i <= start_after:
                continue

        if start_from is not None:
            if not isinstance(start_from, int) or i < start_from:
                continue

        log.info(str(i) + '. ' + step.name)
        if step.run(overwrite, ask_before) == 1:
            log.warning('\n   Process was not complete. You can restart from this point '
                        'using --start-with "' + step.name + '"')
            return 1
        log.info('Done ' + step.name.lower())
        log.info('')
    return 0


if __name__ == '__main__':
    params = parse_args(sys.argv[1:])

    if not isdir(params.out_dir):
        mkdir(params.out_dir)

    with open(params.species_file) as sf:
        species_names = [l.strip() for l in sf.readlines()]
    if isfile(join(params.out_dir, basename(params.species_file))):
        remove(join(params.out_dir, basename(params.species_file)))
    copy(params.species_file, params.out_dir)

    start_after = None
    start_from = params.start_from
    if params.start_from == 'uselog':
        with open(join(params.out_dir, log_file)) as f:
            for l in f.readlines():
                if len(l) > 41 and l[36:40] == 'Done':
                    start_after = l[41:].strip()
                    start_from = None

    set_up_logging(params.debug, params.out_dir)

    log.info('Changing to %s' % params.out_dir)
    chdir(params.out_dir)

    run_workflow(params.out_dir,
                 species_names,
                 params.overwrite or start_after != 0 or start_from != 0,
                 params.blast_dbsize,
                 params.ask_each_step,
                 start_after,
                 start_from,
                 params.threads,
                 params.proxy)