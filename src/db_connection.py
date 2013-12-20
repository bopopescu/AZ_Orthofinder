from os import system
import subprocess
import mysql.connector
from mysql.connector import errorcode
from os.path import join, dirname, realpath

orthomcl_config = join(dirname(realpath(__file__)), 'orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')

import config
import logging
log = logging.getLogger(config.log_fname)


class DbCursor:
    def __init__(self):
        with open(orthomcl_config) as f:
            conf = dict(l.strip().split('=', 1) for l in f.readlines() if l.strip()[0] != '#')
        self.db_login = conf['dbLogin']
        self.db_passw = conf['dbPassword']

        self.cnx = None
        self.cursor = None

    def __connect(self):
        self.cnx = mysql.connector.connect(
            user=self.db_login,
            password=self.db_passw,
            host='127.0.0.1',
            port=3307,
            database='orthomcl',
            buffered=True)

    def __enter__(self):
        connected = False
        while not connected:
            try:
                self.__connect()
                connected = True
            except mysql.connector.errors.InterfaceError:
                #log.info('   MySql server must not be running, trying to start.')
                with open(config.config) as cf:
                    conf = dict(l.strip().split('=', 1) for l
                                in cf.readlines() if l.strip()[0] != '#')

                cmd = 'mysqld_safe --port=%s &' % conf['db_port']
                log.info('   WARNING: Could not connect to MySql server. If it is not running, '
                         'please, start it with mysqld_safe --port=%s &' % conf['db_port'])
                log.info('   Then press any key to proceed.')
                try:
                    raw_input('')
                except KeyboardInterrupt:
                    exit(1)
                #log.info('   ' + cmd)
                #result = subprocess.call(cmd.split())
                #print result
                #if result == 0:
                #    self.__connect()

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