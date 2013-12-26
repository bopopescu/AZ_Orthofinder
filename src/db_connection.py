from os import system
import subprocess
import mysql.connector
from mysql.connector import errorcode
from os.path import join, dirname, realpath

import config
import logging
log = logging.getLogger(config.log_fname)


class DbCursor:
    def __init__(self, step=''):
        with open(config.config_file) as f:
            conf = dict(l.strip().split('=', 1) for l in f.readlines() if l.strip()[0] != '#')
        self.db_login = conf['db_login']
        self.db_passw = conf['db_password']
        self.db_port = conf['db_port']

        self.cnx = None
        self.cursor = None

        self.step = step

    def __connect(self):
        self.cnx = mysql.connector.connect(
            user=self.db_login,
            password=self.db_passw,
            host='127.0.0.1',
            port=self.db_port,
            database='orthomcl',
            buffered=True)

    def __enter__(self):
        connected = False
        while not connected:
            try:
                self.__connect()
                connected = True

            except Exception:
                #log.info('   MySql server must not be running, trying to start.')
                with open(config.config_file) as cf:
                    conf = dict(l.strip().split('=', 1) for l
                                in cf.readlines() if l.strip()[0] != '#')

                cmd = 'mysqld_safe --port=%s &' % conf['db_port']
                log.info('   Could not connect to MySql server. If it is not running, '
                         'please, start it at another terminal with "mysqld_safe --port=%s &"' % conf['db_port'])
                try:
                    raw_input('   After that, press any key to proceed to the next step, or type Ctrl-C to quit. '
                              '(Notice that you can start from this step using --start-from%s): '
                              % (' "' + self.step + '"' if self.step else ''))
                except KeyboardInterrupt:
                    print ''
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