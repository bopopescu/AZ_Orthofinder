from os import system
import subprocess
import mysql_python.connector
from mysql_python.connector import errorcode
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
        self.db_server = conf['db_server']

        self.cnx = None
        self.cursor = None

        self.step = step

    def __connect(self):
        self.cnx = mysql_python.connector.connect(
            user=self.db_login,
            password=self.db_passw,
            host=self.db_server,
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

                cmd1 = 'mysqld --port=%s &' % conf['db_port']
                cmd2 = 'mysqld_safe --port=%s &' % conf['db_port']
                log.info('   Could not connect to MySql server. '
                         'The connection address: %s:%s, login: %s, password: %s, database: %s. '
                         'Please, review the connection options in the confix.txt file in the root directory.' %
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

        self.cursor = self.cnx.cursor()
        return self.cursor

    def __exit__(self, type, err, traceback):
        self.cursor.close()
        self.cnx.close()
        if isinstance(err, mysql_python.connector.Error):
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                log.error('Either incorrect user name or password')
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                log.error('Database orthomcl does not exist')
            else:
                log.error(err)