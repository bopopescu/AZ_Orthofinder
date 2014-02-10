from os import system
import subprocess
import traceback
import sqlite3
import mysql.connector
from mysql.connector import errorcode
from os.path import join, dirname, realpath

import config
import logging
log = logging.getLogger(config.log_fname)

import sys
sys.path = [config.src_dir] + sys.path


class DbCursor:
    def __init__(self, data_fpath='', step=''):
        with open(config.config_file) as f:
            conf = dict(l.strip().lower().split('=', 1) for l
                        in f.readlines() if l.strip()[0] != '#')

        self.db_vendor = conf['db_vendor']

        if self.db_vendor == 'mysql':
            self.db_login = conf['db_login']
            self.db_passw = conf['db_password']
            self.db_port = conf['db_port']
            self.db_server = conf['db_server']
        else:
            self.data_fpath = data_fpath

        self.conn = None
        self.cursor = None

        self.step = step

    def __connect(self):
        if self.db_vendor == 'mysql':
            self.conn = mysql.connector.connect(
                user=self.db_login,
                password=self.db_passw,
                host=self.db_server,
                port=self.db_port,
                database='orthomcl',
                buffered=True)
            self.conn.autocommit = True
        else:
            self.conn = sqlite3.connect(self.data_fpath)
            self.conn.autocommit = True

    def __enter__(self):
        if self.db_vendor == 'mysql':
            connected = False
            while not connected:
                try:
                    self.__connect()
                    connected = True

                except mysql.connector.Error, err:
                    #log.info('   MySql server must not be running, trying to start.')
                    with open(config.config_file) as cf:
                        conf = dict(l.strip().split('=', 1) for l
                                    in cf.readlines() if l.strip()[0] != '#')

                    log.warn('   ' + str(traceback.format_exc(limit=0)))

                    cmd1 = 'mysqld --port=%s &' % conf['db_port']
                    cmd2 = 'mysqld_safe --port=%s &' % conf['db_port']
                    log.info('   Could not connect to MySql server. '
                             'The connection address: %s, port: %s, login: %s, password: %s, database: %s. '
                             'Please, review the connection options in the config.txt file in the root directory.' %
                             (self.db_server, self.db_port, self.db_login, self.db_passw, 'orthomcl'))
                    log.info('   If the server is not running, '
                             'please, start it at another terminal with "%s" or "%s"' % (cmd1, cmd2))
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

                except sqlite3.OperationalError, err:
                    log.exception('SQLite DB error. Please, send the log file to vlad.saveliev@me.com.')
                    return 1

        else:
            self.__connect()

        self.cursor = self.conn.cursor()
        return self.cursor

    def __exit__(self, type, err, traceback):
        self.cursor.close()
        self.conn.close()

        if self.db_vendor == 'mysql':
            if isinstance(err, mysql.connector.Error):
                if err.errno == mysql.connector.errorcode.ER_ACCESS_DENIED_ERROR:
                    log.error('Either incorrect user name or password')
                elif err.errno == mysql.connector.errorcode.ER_BAD_DB_ERROR:
                    log.error('Database orthomcl does not exist')
                else:
                    log.error(err)