from genericpath import isfile, isdir
from itertools import izip, count, ifilterfalse, ifilter
from os import remove
from os.path import basename, exists
from shutil import rmtree
import sqlite3
import subprocess
from db_connection import DbCursor
from mysql.connector import errorcode
import mysql

import config
import logging
log = logging.getLogger(config.log_fname)

orthomcl_config = config.orthomcl_config
orthomcl_bin_dir = config.orthomcl_bin_dir


class Workflow:
    def __init__(self, working_dir, id, cmdline_args):
        self.working_dir = working_dir
        self.id = id
        self.steps = []
        self.cmdline_args = cmdline_args

    def add(self, step):
        self.steps.append(step)

    def extend(self, steps):
        self.steps.extend(steps)

    def run(self, start_after, start_from, overwrite, ask_before):
        if start_from is not None \
                and isinstance(start_from, basestring):
            step_found = False
            for step in self.steps:
                step_found = (start_from.lower() == step.name.lower())
            if not step_found:
                log.error('Step "%s" was not found, maybe incorrect spelling? The list of available steps:' % start_from)
                for i, step in izip(count(1), filter(None, self.steps)):
                    log.error('  %d. %s' % (i, step.name))
                log.error('')

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
            res = step._run(overwrite, ask_before)
            if res != 0:
                if '--start-from' in self.cmdline_args:
                    inx = self.cmdline_args.index('--start-from')
                    del self.cmdline_args[inx]
                    del self.cmdline_args[inx]

                self.cmdline_args.append('--start-from')
                self.cmdline_args.append(str(i))

                log.warn('')
                log.warn('   Process was not complete. You can restart from this point with the following:')
                log.warn('      ' + ' '.join(self.cmdline_args))
                return 1

            log.info('   Done.')
            log.info('')
        return 0


def cmdline(command, parameters=None, stdin=None,
            stdout='pipe', stderr='pipe', env=None,
            ignore_lines_by_pattern=None,
            start_ignoring_from=None):
    parameters = parameters or []

    def callback():
        commandline = ' '.join([command] + map(str, parameters))

        stdin_f = None
        if stdin:
            commandline += ' < ' + stdin
            stdin_f = open(stdin)

        stdout_f = None
        stderr_f = None

        if stdout:
            if stdout in ['pipe', 'log']:
                stdout_f = subprocess.PIPE
            else:
                stdout_f = open(stdout, 'w')
                commandline += ' > ' + stdout

        if stderr:
            if stderr == stdout:
                stderr_f = subprocess.STDOUT
            elif stderr in ['pipe', 'log']:
                stderr_f = subprocess.PIPE
            else:
                stderr_f = open(stderr, 'w')
                commandline += ' 2> ' + stderr

        log.info('   ' + commandline)
        try:
            p = subprocess.Popen(
                [command] + map(str, parameters), env=env,
                stdin=stdin_f, stdout=stdout_f, stderr=stderr_f)

            if stdout_f == subprocess.PIPE:
                for line in iter(p.stdout.readline, ''):

                    if start_ignoring_from:
                        import re
                        a = re.compile(start_ignoring_from)
                        if a.match(line.strip()):
                            break

                    if ignore_lines_by_pattern:
                        import re
                        a = re.compile(ignore_lines_by_pattern)
                        if a.match(line.strip()):
                            continue

                    if stdout == 'pipe':
                        log.info('   ' + line.strip())
                    if stdout == 'log':
                        log.debug('   ' + line.strip())

            if stderr_f == subprocess.PIPE:
                for line in iter(p.stderr.readline, ''):

                    if start_ignoring_from:
                        import re
                        a = re.compile(start_ignoring_from)
                        if a.match(line.strip()):
                            break

                    if ignore_lines_by_pattern:
                        import re
                        a = re.compile(ignore_lines_by_pattern)
                        if a.match(line.strip()):
                            continue
                            
                    if stderr == 'pipe':
                        log.info('   ' + line.strip())
                    if stderr == 'log':
                        log.debug('   ' + line.strip())
            ret_code = p.wait()
            log.debug('      Ret ' + str(ret_code))
            return ret_code

        except KeyboardInterrupt:
            return 1

        except:
            return 1

    return callback


class Step:
    def __init__(self, name, run,
                 req_files=None, prod_files=None,
                 req_tables=None, prod_tables=None):
        self.name = name
        self.run = run
        self.req_files = req_files or []
        self.req_tables = req_tables or []
        self.prod_files = prod_files or []
        self.prod_tables = prod_tables or []


    #def __run(self):
        #if hasattr(self.command, '__call__'):
        #    return self.command(*self.parameters)
        #else:


    def __check_existence(self, overwrite):
        missing_prod_files = list(ifilterfalse(exists, self.prod_files))

        missing_prod_tables = []
        existing_prod_tables = []
        if self.prod_tables:
            with DbCursor() as cursor:
                for table in self.prod_tables:
                    try:
                        query = 'select count(*) from %s;' % table
                        cursor.execute(query)
                        log.debug('   %s exists' % table)
                    except (mysql.connector.Error, sqlite3.OperationalError):
                        #log.debug('   err.errno == errorcode.ER_TABLE_EXISTS_ERROR: ' +
                        #          str(err.errno == errorcode.ER_TABLE_EXISTS_ERROR))
                        #log.debug(err.msg)
                        missing_prod_tables.append(table)
                        log.debug('   %s does not exist' % table)
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
                return False, existing_prod_tables, 0

        return True, existing_prod_tables, 0


    def __check_requirements(self, overwrite, existing_prod_tables):
        missing_req_files = list(ifilterfalse(exists, self.req_files))
        if missing_req_files:
            log.error('   ' + self.name + ' requires files ' +
                      ', '.join(missing_req_files))
            return False, 1

        missing_req_tables = []
        if self.req_tables:
            with DbCursor() as cursor:
                for table in self.req_tables:
                    try:
                        cursor.execute('select count(*) from %s;' % table)
                    except (mysql.connector.Error, sqlite3.OperationalError):
                        missing_req_tables.append(table)
        if missing_req_tables:
            log.error('   ' + self.name + ' requires tables ' +
                      ', '.join(missing_req_tables) + ' installed')
            return False, 1

        # Removing existing data if overwrite
        existing_prod_files = list(ifilter(exists, self.prod_files))
        if overwrite and existing_prod_files:
            log.info('   Overwriting existing ' + ', '.join(existing_prod_files))
            for file in existing_prod_files:
                if isfile(file):
                    remove(file)
                if isdir(file):
                    rmtree(file)

        log.debug('2 existing_prod_tables = ' + str(existing_prod_tables))
        log.debug('2 overwrite = ' + str(overwrite))
        if overwrite and existing_prod_tables:
            with DbCursor() as cursor:
                for table in existing_prod_files:
                    try:
                        log.debug('   drop table %s;' % table)
                        cursor.execute('drop table %s;' % table)
                    except sqlite3.OperationalError, err:
                        log.critical(err)
                    except mysql.connector.Error, err:
                        log.critical(err)
        return True, 0


    def _run(self, overwrite=False, step_by_step=False):
        if step_by_step:
            raw_input('   Proceed?')

        # Checking existence of produced tables and files
        ok, existing_prod_tables, code = self.__check_existence(overwrite)
        log.debug('existing_prod_tables = ' + str(existing_prod_tables))
        if not ok:
            return code

        # Checking requirements
        log.debug('overwrite = ' + str(overwrite))
        ok, code = self.__check_requirements(overwrite, existing_prod_tables)
        if not ok:
            return code

        return self.run()