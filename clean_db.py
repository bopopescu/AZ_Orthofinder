#!/usr/bin/env python

from os.path import join, dirname, realpath
import sqlite3
from src import config
import src.mysql.connector
from src.mysql.connector import errorcode
import sys

from src.db_connection import DbCursor

orthomcl_config = join(dirname(realpath(__file__)), 'src/orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'src/orthomcl_software/bin')


def print_if_main(txt):
    if __name__ == '__main__':
        print txt


def clean_db(suffixes):
    if suffixes:
        if not isinstance(suffixes, (list, tuple)):
            suffixes = [suffixes]
    else:
        suffixes = ['']

    for suffix in suffixes:
        with open(orthomcl_config) as f:
            conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')

        tables = [t.strip() + suffix for t in [
            'BestHit',
            'BestQueryTaxonScore',
            'BestInterTaxonScore',
            'BetterHit',
            'CoOrthNotOrtholog',
            'CoOrthologTaxon',
            'CoOrthologCandidate',
            'CoOrthologAvgScore',
            'CoOrthologTemp',
            'CoOrtholog',
            'InParalog2Way',
            'InParalogAvgScore',
            'InParalog',
            'InParalogTemp',
            'InParalogTaxonAvg',
            'InParalogOrtholog',
            'InplgOrthTaxonAvg',
            'InplgOrthoInplg',
            'InterTaxonMatch',
            'Ortholog',
            'OrthologAvgScore',
            'OrthologTemp',
            'Ortholog2Way',
            'OrthologTaxon',
            'OrthologUniqueId',
            'UniqSimSeqsQueryId',
            'SimilarSequences']]

        with DbCursor(data_fpath=config.sqlite_file) as cursor:
            for table in tables:
                try:
                    query = 'drop table %s;' % table
                    print_if_main(query)
                    cursor.execute(query)
                except src.mysql.connector.Error, err:
                    print_if_main(err.msg)
                    pass
                except sqlite3.OperationalError, err:
                    print_if_main(str(err))
                    pass
            try:
                query = 'drop view %s;' % (conf['interTaxonMatchView'].strip() + suffix)
                print_if_main(query)
                cursor.execute(query)
            except src.mysql.connector.Error, err:
                print_if_main(err.msg)
                pass
            except sqlite3.OperationalError, err:
                print_if_main(str(err))
                pass
    return 0

if __name__ == '__main__':
    clean_db(sys.argv[1:])


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