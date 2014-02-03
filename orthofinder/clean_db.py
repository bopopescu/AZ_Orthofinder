#!/usr/bin/env python

from os.path import join, dirname, realpath
import src.mysql.connector
from src.mysql.connector import errorcode
import sys

from src.db_connection import DbCursor

orthomcl_config = join(dirname(realpath(__file__)), 'src/orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'src/orthomcl_software/bin')


def prt(txt):
    if __name__ == '__main__':
        print txt


def clean_db(suffixes):
    if not isinstance(suffixes, (list, tuple)):
        suffixes = [suffixes]

    for suffix in suffixes:
        with open(orthomcl_config) as f:
            conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')

        tables = [t.strip() + suffix for t in [
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
                    prt(query)
                    cursor.execute(query)
                except src.mysql.connector.Error, err:
                    prt(err.msg)
                    pass
            try:
                query = 'drop view %s;' % (conf['interTaxonMatchView'].strip() + suffix)
                prt(query)
                cursor.execute(query)
            except src.mysql.connector.Error, err:
                prt(err.msg)
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