from genericpath import isfile, isdir
from itertools import izip, count, ifilterfalse, ifilter
from os import remove, getcwd
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
            starting_from_here = i == start_after + 1 if start_after else i == start_from if start_from else None
            res = step._run(i, overwrite, ask_before, starting_from_here)
            if res != 0:
                if '--start-from' in self.cmdline_args:
                    inx = self.cmdline_args.index('--start-from')
                    del self.cmdline_args[inx]
                    del self.cmdline_args[inx]

                self.cmdline_args.append('--start-from')
                self.cmdline_args.append(str(i))

                log.warn('')
                log.warn('   The pipeline is not complete. You can use intermediate results in ' + getcwd() +
                         ', or restart from this point using the --start-from option:')
                log.warn('   ' + ' '.join(self.cmdline_args))
                return 1

            log.info('   Done.')
            log.info('')
        return 0


def cmdline(command, parameters=None, stdin=None,
            stdout='pipe', stderr='pipe', env=None,
            ignore_output_lines_by_pattern=None,
            start_ignoring_from=None):
    parameters = parameters or []

    def callback():
        if isinstance(command, basestring):
            command_list = command.split()
        else:
            command_list = command
        command_list = command_list + map(str, parameters)
        command_str = ' '.join(command_list)

        stdin_f = None
        if stdin:
            command_str += ' < ' + stdin
            stdin_f = open(stdin)

        stdout_f = None
        stderr_f = None

        if stdout:
            if stdout in ['pipe', 'log']:
                stdout_f = subprocess.PIPE
            else:
                stdout_f = open(stdout, 'w')
                command_str += ' > ' + stdout

        if stderr:
            if stderr == stdout:
                stderr_f = subprocess.STDOUT
            elif stderr in ['pipe', 'log']:
                stderr_f = subprocess.PIPE
            else:
                stderr_f = open(stderr, 'w')
                command_str += ' 2> ' + stderr

        log.info('   ' + command_str)
        try:
            p = subprocess.Popen(command_list, env=env,
                stdin=stdin_f, stdout=stdout_f, stderr=stderr_f)

            if stdout_f == subprocess.PIPE:
                for line in iter(p.stdout.readline, ''):

                    if start_ignoring_from:
                        import re
                        a = re.compile(start_ignoring_from)
                        if a.match(line.strip()):
                            break

                    if ignore_output_lines_by_pattern:
                        import re
                        a = re.compile(ignore_output_lines_by_pattern)
                        if a.match(line.strip()):
                            continue

                    if stdout == 'pipe':
                        log.info('   ' + line.strip())
                    if stdout == 'log':
                        log.debug('   ' + line.strip())

            stderr_output = []
            if p.stderr:
                for line in iter(p.stderr.readline, ''):
                    if start_ignoring_from:
                        import re
                        a = re.compile(start_ignoring_from)
                        if a.match(line.strip()):
                            break

                    if ignore_output_lines_by_pattern:
                        import re
                        a = re.compile(ignore_output_lines_by_pattern)
                        if a.match(line.strip()):
                            continue

                    if stderr == 'pipe':
                        log.info('   ' + line.strip())
                    else:
                        stderr_output.append(line)
                    if stderr == 'log':
                        log.debug('   ' + line.strip())

            ret_code = p.wait()
            if ret_code != 0:
                log.error('')
                log.error('Command returned ' + str(ret_code))
                for line in stderr_output:
                    log.error('   ' + line.strip())
                    log.error('')
            return ret_code

        except KeyboardInterrupt:
            return 1

        except OSError, e:
            log.error('')
            log.error('   OS Error when executing: ' + command_str)
            log.error('   ' + e.strerror)
            if e.filename:
                log.error('   For ' + e.filename)
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

        with open(config.config_file) as f:
            conf = dict(l.strip().lower().split('=', 1) for l
                        in f.readlines() if l.strip() and l.strip()[0] != '#')
        db_vendor = conf['db_vendor']

        missing_prod_tables = []
        existing_prod_tables = []
        if self.prod_tables:
            with DbCursor() as cursor:
                for table in self.prod_tables:
                    #if db_vendor == 'sqlite':
                    #    try:
                    #        cursor.execute("SELECT name FROM sqlite_master "
                    #                       "WHERE type='table' AND name='%s';" % table)
                    #        if cursor.fetchone():
                    #            existing_prod_tables.append(table)
                    #            log.debug('   %s exists' % table)
                    #        else:
                    #            missing_prod_tables.append(table)
                    #            log.debug('   %s does not exist' % table)
                    #            continue
                    #    except sqlite3.OperationalError:
                    #        missing_prod_tables.append(table)
                    #        log.debug('   %s does not exist' % table)
                    #    else:
                    #        existing_prod_tables.append(table)
                    #        log.debug('   %s exists' % table)
                    #else:
                    try:
                        q = 'SELECT 1 FROM %s LIMIT 1;' % table
                        log.debug(q)
                        cursor.execute(q)
                        res = cursor.fetchall()
                        log.debug(res)
                    except (mysql.connector.Error, sqlite3.OperationalError):
                        missing_prod_tables.append(table)
                        log.debug('   %s does not exist' % table)
                    else:
                        existing_prod_tables.append(table)
                        log.debug('   %s exists' % table)
            log.debug('')

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

        with open(config.config_file) as f:
            conf = dict(l.strip().lower().split('=', 1) for l
                        in f.readlines() if l.strip() and l.strip()[0] != '#')
        db_vendor = conf['db_vendor']

        missing_req_tables = []
        if self.req_tables:
            with DbCursor() as cursor:
                for table in self.req_tables:
                    #if db_vendor == 'sqlite':
                    #    try:
                    #        q = "SELECT name FROM sqlite_master " + \
                    #            "WHERE type='table' AND name='%s';" % table
                    #        log.debug(q)
                    #        cursor.execute(q)
                    #        res = cursor.fetchall()
                    #        log.debug(res)
                    #        if not res:
                    #            missing_req_tables.append(table)
                    #            continue
                    #    except sqlite3.OperationalError:
                    #        missing_req_tables.append(table)
                    #else:
                    try:
                        q = 'SELECT 1 FROM %s LIMIT 1;' % table
                        log.debug(q)
                        cursor.execute(q)
                        res = cursor.fetchall()
                        log.debug(res)
                    except (mysql.connector.Error, sqlite3.OperationalError):
                        log.exception('aaa')
                        missing_req_tables.append(table)

        if missing_req_tables:
            log.error('   ' + self.name + ' requires tables ' +
                      ', '.join(missing_req_tables))
            return False, 1

        # Removing existing data if overwrite
        existing_prod_files = list(ifilter(exists, self.prod_files))
        if overwrite and existing_prod_files:
            log.info('   overwriting ' + ', '.join(existing_prod_files))
            for file in existing_prod_files:
                if isfile(file):
                    remove(file)
                if isdir(file):
                    rmtree(file)

        if overwrite and existing_prod_tables:
            with DbCursor() as cursor:
                for table in existing_prod_tables:
                    try:
                        log.debug('   drop table %s;' % table)
                        cursor.execute('drop table %s;' % table)
                    except sqlite3.OperationalError, err:
                        log.critical(err)
                    except mysql.connector.Error, err:
                        log.critical(err)
        return True, 0


    def _run(self, step_number, overwrite=False, step_by_step=False, starting_from_here=False):
        if step_by_step:
            raw_input('   Proceed?')

        # Checking existence of produced tables and files
        ok, existing_prod_tables, code = self.__check_existence(overwrite)
        if not ok:
            return code

        # Checking requirements
        ok, code = self.__check_requirements(overwrite, existing_prod_tables)
        if not ok:
            return code

        return self.run(starting_from_here=starting_from_here)