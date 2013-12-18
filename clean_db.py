from os.path import join, dirname, realpath
from mysql.connector import errorcode
import mysql.connector
import sys

from src.db_connection import DbCursor

orthomcl_config = join(dirname(realpath(__file__)), 'src/orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'src/orthomcl_software/bin')


def clean_db(workflow_id):
    with open(orthomcl_config) as f:
        conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')

    tables = [t.strip() + '_' + workflow_id for t in [
        conf['orthologTable'],
        conf['inParalogTable'],
        conf['coOrthologTable'],
        conf['similarSequencesTable'],
        'BestInterTaxonScore',
        'CoOrthNotOrtholog',
        'CoOrthologTaxon',
        'CoOrthologCandidate',
        'CoOrthologAvgScore',
        'CoOrthologTemp',
        'BetterHit',
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
        'BestHit',
        'BestQueryTaxonScore']]

    with DbCursor() as cursor:
        for table in tables:
            try:
                query = 'drop table %s;' % table
                print query
                cursor.execute(query)
            except mysql.connector.Error, err:
                print err.msg
                pass
        try:
            query = 'drop view %s;' % (conf['interTaxonMatchView'].strip() + '_' + workflow_id)
            print query
            cursor.execute(query)
        except mysql.connector.Error, err:
            print err.msg
            pass


if __name__ == '__main__':
    clean_db(sys.argv[1])


#    #cmd = ' '.join([
#    #    join(orthomcl_bin_dirpath, 'orthomclPairs'),
#    #    orthomcl_config_fpath,
#    #    'orthomclpairs.log',
#    #    'cleanup=all',
#    #    ('suffix=_' + workflow_id) if workflow_id else '',
#    #    ])
#    #print '   ' + cmd
#    #system(cmd)
#
#    cmd = ' '.join([
#        join(orthomcl_bin_dir, 'orthomclDropSchema'),
#        orthomcl_config,
#        'drop_log.sql',
#        ('_' + workflow_id) if workflow_id else '',
#        ])
#    log.info('   ' + cmd)
#    system(cmd)
#    log.info()