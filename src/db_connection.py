from os.path import join, dirname, realpath
import mysql.connector
from mysql.connector import errorcode

orthomcl_config = join(dirname(realpath(__file__)), 'orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')

import logging
log = logging.getLogger('orthofinder')


class DbCursor:
    def __init__(self):
        with open(orthomcl_config) as f:
            conf = dict(l.split('=') for l in f.readlines() if l[0] != '#')
        self.db_login = conf['dbLogin']
        self.db_passw = conf['dbPassword']
        self.cnx = None
        self.cursor = None

    def __enter__(self):
        self.cnx = mysql.connector.connect(
            user=self.db_login,
            password=self.db_passw,
            host='127.0.0.1',
            port=3307,
            database='orthomcl',
            buffered=True)

        self.cursor = self.cnx.cursor()
        return self.cursor

    def __exit__(self, type, err, traceback):
        self.cursor.close()
        self.cnx.close()
        if isinstance(err, mysql.connector.Error):
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                log.error('Either incorrect user name or password')
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                log.error('Database orthomcl does not exist')
            else:
                log.error(err)