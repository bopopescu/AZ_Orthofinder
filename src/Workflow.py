from genericpath import isfile, isdir
from itertools import izip, count, ifilterfalse, ifilter
from os import remove
from os.path import basename, exists
from shutil import rmtree
from subprocess import call
import subprocess
from mysql.connector import errorcode
from db_connection import DbCursor, mysql

import config
import logging
log = logging.getLogger(config.log_fname)

orthomcl_config = config.orthomcl_config
orthomcl_bin_dir = config.orthomcl_bin_dir


class Workflow:
    def __init__(self, working_dir, id):
        self.working_dir = working_dir
        self.id = id
        self.steps = []

    def add(self, step):
        self.steps.append(step)

    def extend(self, steps):
        self.steps.extend(steps)

    def run(self, start_after, start_from, overwrite, ask_before):
        for i, step in izip(count(1), filter(None, self.steps)):
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
            if step.run(overwrite, ask_before) != 0:
                log.info('')
                log.warning('   Process was not complete. You can restart from this point '
                            'using --start-from "' + step.name + '"')
                return 1
            log.info('   Done ' + step.name.lower())
            log.info('')

        return 0


class Step:
    def __init__(self, name, cmd,
                 req_files=None, prod_files=None,
                 req_tables=None, prod_tables=None,
                 parameters=None, stdin=None,
                 stdout='pipe', stderr='pipe'):
        self.name = name
        self.command = cmd
        self.req_files = req_files or []
        self.req_tables = req_tables or []
        self.prod_files = prod_files or []
        self.prod_tables = prod_tables or []
        self.parameters = parameters or []
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

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
                        query = '   select count(*) from %s;' % table
                        cursor.execute(query)
                    except mysql.connector.Error, err:
                        #log.debug('   err.errno == errorcode.ER_TABLE_EXISTS_ERROR: ' +
                        #          str(err.errno == errorcode.ER_TABLE_EXISTS_ERROR))
                        #log.debug(err.msg)
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

            stdin_f = None
            if self.stdin:
                commandline += ' < ' + self.stdin
                stdin_f = open(self.stdin)

            stdout_f = subprocess.PIPE
            stderr_f = subprocess.PIPE
            if self.stdout and self.stdout not in ['log', 'pipe']:
                stdout_f = open(self.stdout, 'w')
                commandline += ' > ' + self.stdout
            else:
                stderr_f = subprocess.STDOUT

            log.info('   ' + commandline)
            try:
                p = subprocess.Popen([self.command] + map(str, self.parameters),
                                     stdin=stdin_f, stdout=stdout_f, stderr=stderr_f)
                if stdout_f == subprocess.PIPE:
                    for line in iter(p.stdout.readline, ''):
                        if self.stdout == 'pipe':
                            log.info('   ' + line.strip())
                        if self.stdout == 'log':
                            log.debug('   ' + line.strip())
                if stderr_f == subprocess.PIPE:
                    for line in iter(p.stderr.readline, ''):
                        if self.stderr == 'pipe':
                            log.info('   ' + line.strip())
                        if self.stderr == 'log':
                            log.debug('   ' + line.strip())
                ret_code = p.wait()
                log.debug('      Ret ' + str(ret_code))
                return ret_code

            except KeyboardInterrupt:
                return 1

            except:
                return 1
