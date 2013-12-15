from genericpath import isfile, isdir, exists
from itertools import ifilterfalse, ifilter
from os import chdir, system
from os.path import join, realpath, dirname, relpath, basename, normpath
from subprocess import call
from sys import stderr
from Bio.Blast.Applications import NcbiblastpCommandline
from fetch_proteomes import fetch_proteomes
import mysql.connector
from mysql.connector import errorcode

file_dirpath = dirname(realpath(__file__))


class DbCursor:
    def __init__(self, orthomcl_config):
        with open(orthomcl_config) as f:
            params = dict(l.split('=') for l in f.readlines() if l[0] != '#')
        self.db_login = params['dbLogin']
        self.db_passw = params['dbPassword']

    def __enter__(self, orthomcl_config):
        self.cnx = mysql.connector.connect(
            user=self.db_login,
            password=self.db_passw,
            host='localhost',
            database='orthomcl')
        return self.cnx.cursor()

    def __exit__(self, type, err, traceback):
        self.cnx.close()
        if isinstance(err, mysql.connector.Error):
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print('Either incorrect user name or password')
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print('Database orthomcl does not exist')
            else:
                print(err)

class Step:
    def __init__(self, name, cmd,
                 req_files=None,
                 req_tables=None,
                 prod_files=None,
                 prod_tables=None,
                 parameters=None,
                 orthomcl_config=None):
        self.name = name
        self.command = cmd
        self.req_files = req_files or []
        self.req_tables = req_tables or []
        self.prod_files = prod_files or []
        self.prod_tables = prod_tables or []
        self.parameters = parameters or []
        self.orthomcl_config = orthomcl_config

    def program(self):
        return basename(self.command)

    def run(self, reuse=True):
        assert self.command

        missing_req_files = list(ifilterfalse(exists, self.req_files))
        if missing_req_files:
            for file in missing_req_files:
                print >> stderr, '   ' + self.name + ' requires ' + file
            return 1

        missing_prod_files = list(ifilterfalse(isfile, self.prod_files))
        if reuse and self.prod_files and not missing_prod_files:
            print '   Skipping: files ' + ' '.join(self.prod_files) + ' produces already exist'
            return 0

        #if self.req_tables:
        #    with DbCursor(self.orthomcl_config) as cursor:
        #        for table in self.req_tables:
        #            try:
        #                cursor.execute('truncate table %s;' % table)
        #            except mysql.connector.Error, err:
        #                print 'Failed deleting data from %s: %s' % (table, err)
        #
        #if self.prod_tables:
        #    with DbCursor(self.orthomcl_config) as cursor:
        #        for table in self.prod_tables:
        #            try:
        #                cursor.execute('truncate table %s;' % table)
        #            except mysql.connector.Error, err:
        #                print 'Failed deleting data from %s: %s' % (table, err)

        raw_input('   Proceed?')
        print '   ' + ' '.join([self.program()] + map(str, self.parameters))
        return system(' '.join([self.command] + map(str, self.parameters)))

def run_workflow(
        working_dirpath,
        orthomcl_bin_dirpath,
        orthomcl_config_fpath,
        workflow_id='',
        reuse=True,
        inflation=1.5,
        start_id_with=1000):

    orthomcl_bin = orthomcl_bin_dirpath
    orthomcl_config = orthomcl_config_fpath
    working_dirpath = normpath(working_dirpath)
    sql_log = 'log.sql'

    proteomes = 'proteomes'
    good_proteins = 'good_proteins.fasta'
    poor_proteins = 'poor_proteins.fasta'
    blast_out = 'blasted.tsv'
    similar_sequences = 'similar_sequences.txt'
    pairs_log = 'orthomclpairs.log'
    groups_file = 'groups.txt'
    singletons_file = 'singletons.txt'

    with open(orthomcl_config) as f:
        params = dict(l.split('=') for l in f.readlines() if l[0] != '#')
        ortholog_table = params['orthologTable']
        in_paralog_table = params['inParalogTable']
        coortholog_table = params['coOrthologTable']
        similar_sequeces_table = params['similarSequencesTable']
        best_hit_table = 'BestHit'
        best_hit_taxon_score_table = 'BestQueryTaxonScore'

    print 'Working at', working_dirpath
    chdir(working_dirpath)
    steps = []

    if not reuse:
        steps.append(Step(
            'Installing schema',
            cmd=join(orthomcl_bin, 'orthomclInstallSchema'),
            req_files=[orthomcl_config],
            prod_files=[],  # drops and creates SimilarSequences, InParalog,
                            # Ortholog, CoOrtholog, InterTaxonMatch view
            parameters=[
                orthomcl_config,
                sql_log,
                '_' + workflow_id,
                ]))

    steps.extend([
        Step('Filtering fasta',
             cmd=join(orthomcl_bin, 'orthomclFilterFasta'),
             req_files=[proteomes],
             prod_files=[good_proteins, poor_proteins],
             parameters=[
                proteomes,
                10,
                20,
                good_proteins,
                poor_proteins
                ]),

        Step('Blasting all vs. all',
             cmd='blastp',
             req_files=[proteomes, good_proteins],
             prod_files=[blast_out],
             parameters=[
                '-query', good_proteins,
                '-subject', good_proteins,
                '-out', blast_out,
                '-outfmt', 6,  # tabular
                '-evalue', 1e-5,
                '-num_descriptions', 10000,  # don't care value
                '-num_alignments', 10000,  # don't care value
                ]),

        Step('Parsing blast results',
             cmd=join(orthomcl_bin, 'orthomclBlastParser'),
             req_files=[proteomes, blast_out],
             prod_files=[similar_sequences],
             parameters=[
                blast_out,
                proteomes,
                '>>', similar_sequences,
                ]),

        Step('Loading blast results into the database',
             cmd=join(orthomcl_bin, 'orthomclLoadBlast'),
             req_files=[orthomcl_config, similar_sequences],  # and initialized database
             prod_files=[],  # loads blast results into the db
             parameters=[
                orthomcl_config,
                similar_sequences,
                ('_' + workflow_id) if workflow_id else '',
                ]),

        Step('Finding pairs',
             cmd=join(orthomcl_bin, 'orthomclPairs'),
             req_files=[orthomcl_config],  # and initialized database
             prod_files=[],  # populates InParalog, Ortholog, CoOrtholog
             parameters=[
                orthomcl_config,
                pairs_log,
                'cleanup=no',
                'startAfter=useLog' if reuse else '',
                ('suffix=_' + workflow_id) if workflow_id else '',
                ]),

        Step('Dump pairs files',
             cmd=join(orthomcl_bin, 'orthomclDumpPairsFiles'),
             req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
             prod_files=[
                'mclInput',
                'pairs/potentialOrthologs.txt',
                'pairs/potentialInparalogs.txt',
                'pairs/potentialCoorthologs.txt',],
             parameters=[
                orthomcl_config,
                working_dirpath,  # no effect, just to use next positional argument
                ('_' + workflow_id) if workflow_id else '']),

        Step('MCL',
             cmd='mcl',
             req_files=['mclInput'],
             prod_files=['mclOutput'],
             parameters=[
                'mclInput',
                '--abc',
                '-I', str(inflation),
                '-o', 'mclOutput',
                ]),

        Step('MCL groups to files',
             cmd=join(orthomcl_bin, 'orthomclMclToGroups'),
             req_files=['mclOutput'],
             prod_files=[groups_file],
             parameters=[
                workflow_id + '_',
                str(start_id_with),
                '<', 'mclOutput',
                '>', groups_file]),

        Step('MCL singletones to files',
             cmd=join(orthomcl_bin, 'orthomclSingletons'),
             req_files=[good_proteins, groups_file],
             prod_files=[singletons_file],
             parameters=[
                good_proteins,
                groups_file,
                '>>', singletons_file]),
    ])
    for i, step in enumerate(steps):
        print str(i + 1) + '.', step.name
        if step.run(reuse) == 1:
            return 1
        print
    return 0

def cleanup(
        working_dirpath,
        orthomcl_bin_dirpath,
        orthomcl_config_fpath,
        workflow_id=None):

    chdir(working_dirpath)
    print 'Cleaning up:'

    cmd = ' '.join([
        join(orthomcl_bin_dirpath, 'orthomclPairs'),
        orthomcl_config_fpath,
        'orthomclpairs.log',
        'cleanup=all',
        ('suffix=_' + workflow_id) if workflow_id else '',
        ])
    print '   ' + cmd
    system(cmd)

    cmd = ' '.join([
        join(orthomcl_bin_dirpath, 'orthomclDropSchema'),
        orthomcl_config_fpath,
        'drop.sql',
        ('_' + workflow_id) if workflow_id else '',
        ])
    print '   ' + cmd
    system(cmd)

#def _clean_db(self, orthomcl_config):
#    with open(orthomcl_config) as f:
#        params = dict(l.split('=') for l in f.readlines() if l[0] != '#')
#
#    db_login = params['dbLogin']
#    db_passw = params['dbPassword']
#    tables = ['BestHit',
#              'BestQueryTaxonScore',
#              params['similarSequencesTable']]
#    views = ['intertaxonmatch']
#
#    import mysql.connector
#    from mysql.connector import errorcode
#    try:
#        cnx = mysql.connector.connect(
#            user=db_login,
#            password=db_passw,
#            host='localhost',
#            database='orthomcl')
#    except mysql.connector.Error as err:
#        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
#            print('Either incorrect user name or password')
#        elif err.errno == errorcode.ER_BAD_DB_ERROR:
#            print('Database orthomcl does not exist')
#        else:
#            print(err)
#    else:
#        cursor = cnx.cursor()
#        for table in tables:
#            try:
#                cursor.execute('truncate table %s;' % table)
#            except mysql.connector.Error, err:
#                print 'Failed deleting data from %s: %s' % (table, err)
#        cnx.close()

    #def __rp(self, path):
    #    return relpath(path, self.species_dirpath)


if __name__ == '__main__':
    orthomcl_bin_dirpath = join(file_dirpath, 'orthomcl_software/bin')
    orthomcl_config_fpath = join(file_dirpath, 'orthomcl.config')

    #species_dirpath = fetch_proteomes('Escherichia coli K-12')
    working_dirpath = join(file_dirpath,
                           '../data/Escherichia_coli_K-12_small_2')
    cleanup(working_dirpath,
            orthomcl_bin_dirpath,
            orthomcl_config_fpath,
            workflow_id='ecoli_0')
    print

    run_workflow(
        working_dirpath,
        orthomcl_bin_dirpath,
        orthomcl_config_fpath,
        workflow_id='ecoli_0',
        reuse=False)
