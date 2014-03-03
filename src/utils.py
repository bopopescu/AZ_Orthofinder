from dircache import listdir
from os import environ, access, X_OK, remove, getcwd, chdir, mkdir
import os
from os.path import basename, isdir, split, join, pathsep, isfile, dirname, exists, devnull, realpath
from random import randint
from shutil import copy
from subprocess import call, Popen, PIPE
import re
import tarfile
import shutil
import config
from config import config_file, orthomcl_config_fname, log_fname, mcl_dir, \
    mysql_cnf, mysql_linux_tar, mysql_osx_tar, mysql_extracted_dir, src_dir
import logging
from Workflow import cmdline

log = logging.getLogger(log_fname)


def make_workflow_id__from_species_names(species_names=None):
    if not species_names:
        return str(randint(1000, 9999))

    sn_words = ' '.join(species_names).split()
    if len(sn_words) == 1:
        return sn_words[0].replace('-', '')[:4]
    if len(sn_words) >= 2:
        return ''.join(w[0] for w in sn_words)


def make_workflow_id(working_dir=None):
    if not working_dir:
        return str(randint(1000, 9999))

    wid = (basename(working_dir) or basename(dirname(working_dir)))
    wid = re.sub('[^0-9a-zA-Z_-]+', '_', wid)
    return wid


def interrupt(msg, code=1):
    log.error(msg)
    exit(code)


def register_ctrl_c():
    from signal import signal, SIGINT
    signal(SIGINT, lambda s, f: interrupt('', 0))


def check_and_install_tools(debug, log_path):
    check_installed_tools(['blastp'], only_warn=True)

    check_install_mcl(debug, log_path=log_path, only_warn=False)

    #prepare_mysql_config()

    check_perl_modules(debug, only_warn=True)


def which(program):
    def is_exe(fpath_):
        return isfile(fpath_) and access(fpath_, X_OK)

    fpath, fname = split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in environ["PATH"].split(pathsep):
            path = path.strip('"')
            exe_file = join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def check_mysql(only_warn=False):
    if which('mysqld'):
        return 'mysqld'
    else:
        if only_warn:
            log.warn('WARNING: MySQL is not installed. It is required for some steps.')
        else:
            log.error('ERROR: MySQL is not installed.')
            exit(1)


def check_install_mysql(only_warn=False, internet_ok=False):
    if which('mysqld'):
        return 'mysqld'
    else:
        if only_warn:
            log.warn('WARNING: MySQL is not installed. It is required for some steps.')
            intenet_ok = test_internet_conn(url='http://dev.mysql.com/get/Downloads/MySQL-5.6/mysql-5.6.15-linux-glibc2.5-x86_64.tar.gz')
            try:
                if intenet_ok:
                    raw_input('Install mysql manually? It will be downloaded from the Internet, '
                              'about 230MB. '
                              'Press any key for "yes", Ctrl-C for "no":')
                else:
                    log.error('We can install mysql for you, but Internet connection seems to be off, '
                              'which is required to download mysql. '
                              'Please, retry with a working Internet connection, '
                              'or install mysql manually.')
                    exit(1)
            except KeyboardInterrupt:
                print ''
                exit(1)

            import urllib
            import tarfile
            link = 'http://dev.mysql.com/get/Downloads/MySQL-5.6/mysql-5.6.15-linux-glibc2.5-x86_64.tar.gz'
            handle = urllib.urlopen(link)
            with open('mysql.tar.gz', 'wb') as out:
                while True:
                    data = handle.read(1024)
                    if len(data) == 0: break
                    out.write(data)

            with tarfile.open('mysql.tar.gz') as tar:
                tar.extractall()

        else:
            log.error('ERROR: MySQL is not installed.')


def check_installed_tools(tools, only_warn=False):
    ok = True
    for tool in tools:
        if not which(tool):
            ok = False
            if only_warn:
                log.warn('WARNING: "' + tool + '" might be not installed. '
                         'It is required to peform some following steps.')
                log.warn('')
            else:
                log.error('ERROR: "' + tool + '" might be not installed. It is required to proceed.')

    if not ok and not only_warn:
        exit(1)
    return ok


def check_install_mcl(debug, log_path=log_fname, only_warn=False):
    if which('mcl'):
        return 'mcl', 0

    mcl_bin_path = join(mcl_dir, 'bin', 'bin', 'mcl')
    if exists(mcl_bin_path):
        return mcl_bin_path, 0
    else:
        if exists(mcl_dir):
            shutil.rmtree(mcl_dir)

    tar_gz = join(src_dir, 'mcl.tar.gz')
    if not isfile(join(src_dir, 'mcl.tar.gz')):
        if not only_warn:
            log.error('Error: no file ' % tar_gz)
            return None, 1
        else:
            log.warning('Warning: no file ' % tar_gz)

    with tarfile.TarFile.open(join(src_dir, 'mcl.tar.gz'), 'r:gz') as mcl_tar_gz:
        mcl_tar_gz.extractall(src_dir)

    for fname in listdir(src_dir):
        if fname.startswith('mcl'):
            shutil.move(join(src_dir, fname), mcl_dir)
            break

    log.debug(log_path)
    log.debug(getcwd())

    def exec_cmdline(command):
        log.info('   ' + command)
        if debug:
            res = cmdline(command.split())()
        else:
            res = cmdline(command.split(), stdout=None, stderr=None)()

        if res != 0:
            log.debug('Running ' + command)
            if only_warn:
                log.warning('WARNING: Cannot find or install mcl. '
                            'It required to perform some steps. '
                            'Try to install it manually: http://micans.org/mcl/src')
            else:
                log.error('ERROR: Cannot find or install mcl. '
                          'Try to install it manually: http://micans.org/mcl/src')
            return None, res

    log.info('Compiling MCL...')
    cur_dir = getcwd()
    log.info('Changing to ' + cur_dir)
    chdir(mcl_dir)
    mcl_path = join(mcl_dir, 'bin')

    exec_cmdline('./configure' + ' -q --prefix=' + mcl_path)
    exec_cmdline('make')
    exec_cmdline('make check')
    exec_cmdline('make install')
    chdir(cur_dir)

    return mcl_bin_path, 0


def check_perl_modules(debug, only_warn=False):
    def check(env=None):
        dbimysql = call(
            'perl -MDBD::mysql -e 1'.split(),
             stderr=open(devnull, 'w'),
             stdout=open(devnull, 'w'),
             env=env) == 0
        if not dbimysql:
            log.debug('DBD::mysql is not installed.')

        dbd = call(
            'perl -MDBI -e 1'.split(),
             stderr=open(devnull, 'w'),
             stdout=open(devnull, 'w'),
             env=env) == 0
        if not dbd:
            log.debug('DBI is not installed.')

        return dbimysql and dbd

    if not check():  # (env):
        if only_warn:
            log.warning('WARNING: Could not find or install Perl modules.')
        else:
            log.error('ERROR: Could not find or install Perl modules. ')

        print '''Pleasy, try the following manually:
    $ perl -MCPAN -e shell
    cpan> install Data::Dumper
    cpan> install DBI
    cpan> force install DBD::mysql
    '''
        if only_warn:
            return 0
        else:
            exit(4)


def check_perl_modules_mysql(debug, only_warn=False):
    def check(env=None):
        dbimysql = call(
            'perl -MDBD::mysql -e 1'.split(),
             stderr=open(devnull, 'w'),
             stdout=open(devnull, 'w'),
             env=env) == 0
        if not dbimysql:
            log.debug('DBD::mysql is not installed.')

        dbd = call(
            'perl -MDBI -e 1'.split(),
             stderr=open(devnull, 'w'),
             stdout=open(devnull, 'w'),
             env=env) == 0
        if not dbd:
            log.debug('DBI is not installed.')

        return dbimysql and dbd

    # tool_patb = join(dirname(realpath(__file__)), '../')
    # assert isdir(tool_patb)
    # perl_modules_path = join(tool_patb, 'perl_modules')
    # lib_path = join(perl_modules_path, 'lib')
    # man_path = join(perl_modules_path, 'man')
    # cpan_path = join(tool_patb, '.cpan')
    # cpan_cpan_path = join(cpan_path, 'CPAN')    #    subprocess.call("perl -wle'print for grep /src/perl_modules/, @INC'",
    #
    # if not isdir(perl_modules_path):
    #     mkdir(perl_modules_path)
    #     mkdir(lib_path)
    #     mkdir(man_path)
    #     mkdir(join(man_path, 'man1'))
    #     mkdir(join(man_path, 'man3'))
    # if not isdir(cpan_path): mkdir(cpan_path)
    # if not isdir(cpan_cpan_path): mkdir(cpan_cpan_path)
    #
    # with open(join(cpan_cpan_path, 'MyConfig.pm'), 'w') as cpan_f, \
    #      open(join(cpan_cpan_path, 'MyConfig_template.pm')) as cpan_tmpl_f:
    #     cpan_f.write(cpan_tmpl_f.read().replace('{PATH}', tool_patb))
    #
    # if 'PERL5LIB' not in environ:
    #     environ['PERL5LIB'] = lib_path
    # else:
    #     environ['PERL5LIB'] = environ['PERL5LIB'] + ':' + lib_path
    #
    # if 'MANPATH' not in environ:
    #     environ['MANPATH'] = lib_path
    # else:
    #     environ['MANPATH'] = environ['MANPATH'] + ':' + man_path
    #
    # env = {'PERL5LIB': lib_path, 'MANPATH': man_path}
    # if check(environ):
    #     return

        #if debug:
        #    subprocess.call("perl -wle'print for grep /src/perl_modules/, @INC'",
        #                    shell=True, env=environ,
        #                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #call('perl -MCPAN -we shell', shell=True,
        #                stdin=open(), stdout=)

        #log.info('Installing Data::Dumper...')
        #cmd = ['cpan', '-j', join(cpan_cpan_path, 'MyConfig.pm'), 'Data::Dumper']
        #log.debug(' '.join(cmd))
        #call(cmd, env=env,
        #     stdout=(open(devnull, 'w') if not debug else PIPE),
        #     stderr=(open(log_file_path, 'w') if not debug else PIPE))

        #log.info('Installing DBI...')
        #cmd = ['cpan', '-j', join(cpan_cpan_path, 'MyConfig.pm'), 'DBI']
        #log.debug(' '.join(cmd))
        #call(cmd, env=env,
        #     stdout=(open(devnull, 'w') if not debug else PIPE),
        #     stderr=(open(log_file_path, 'w') if not debug else PIPE))

    # log.info('Installing DBD::mysql...')
    # # log.debug(' '.join(cmd))
    # cmdline('cpan',
    #         ['-f',
    #          '-rj',
    #          join(cpan_cpan_path, 'MyConfig.pm'),
    #          'DBD::mysql'],
    #         env=env,
    #         stdout='pipe',
    #         stderr='pipe')()

        #p = Popen('perl -MCPAN -we shell'.split(),
        #          stdout=(open(devnull, 'w') if not debug else PIPE),
        #          stdin=PIPE,
        #          stderr=(open(log_file_path, 'w') if not debug else PIPE))
        #
        #grep_stdout = p.communicate(input='one\ntwo\nthree\nfour\nfive\nsix\n')[0]
        #print(grep_stdout)


        #log.info('Installing DBD::mysql perl module...')
        #if which('yum'):
        #    res = call('yum install "perl(DBD::mysql)"'.split(),
        #                stderr=open(devnull, 'w'),
        #                stdout=open(devnull, 'w'))
        #elif which('apt-get'):
        #    res = call('sudo apt-get install libdbd-mysql-perl"'.split(),
        #                stderr=open(devnull, 'w'),
        #                stdout=open(devnull, 'w'))
        #else:
        #    log.error('Cannot install DBD::mysql for your system.')

    if not check():  # (env):
        if only_warn:
            log.warning('WARNING: Could not find or install Perl modules.')
        else:
            log.error('ERROR: Could not find or install Perl modules. ')

        print '''Pleasy, try the following manually:
    $ perl -MCPAN -e shell
    cpan> o conf makepl_arg "mysql_config=%s"
    cpan> install Data::Dumper
    cpan> install DBI
    cpan> force install DBD::mysql
    ''' % mysql_cnf
        if only_warn:
            return 0
        else:
            exit(4)


def read_list(file, where_to_save=None):
    if not file:
        return None
    with open(file) as f:
        results = [l.strip() for l in f.readlines() if l.strip()]
    if where_to_save:
        if isfile(join(where_to_save, basename(file))):
            remove(join(where_to_save, basename(file)))
        copy(file, where_to_save)
    return results


def test_entrez_conn():
    return test_internet_conn('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')


def test_blast_conn():
    return test_internet_conn('http://blast.ncbi.nlm.nih.gov/Blast.cgi')


def test_ftp_conn():
    return test_internet_conn('ftp.ncbi.nih.gov')


def test_internet_conn(url):
    import urllib2
    try:
        response = urllib2.urlopen(url, timeout=1)
    except urllib2.URLError as err:
        log.debug('test_internet_conn fail. Url=' + url + ', error=' + str(err))
        return False
    else:
        return True


def set_up_config(working_dir):
    with open(config_file) as cf:
        conf = dict(
            l.strip().split('=', 1) for l
            in cf.readlines() if l.strip() and l.strip()[0] != '#')
        log.debug('Read conf: ' + str(conf))

    with open(join(src_dir, orthomcl_config_fname)) as ocf:
        omcl_conf = dict(l.strip().split('=', 1) for l
                         in ocf.readlines() if l.strip() and l.strip()[0] != '#')

    memory = int(conf['memory'])

    if conf['db_vendor'] == 'sqlite':
        db_file = join(working_dir, config.sqlite_file)

        omcl_conf['dbVendor'] = 'sqlite'
        omcl_conf['dbConnectString'] = \
            'DBI:SQLite:%s' % db_file

    else:
        omcl_conf['dbVendor'] = 'mysql'
        omcl_conf['dbConnectString'] = \
            'dbi:mysql:database=orthomcl;' \
            'host=%s;' \
            'port=%s;' \
            'myisam_sort_buffer_size=%dG;' \
            'read_buffer_size=%dG;' \
            'innodb_buffer_pool_size=%dG;' \
            'mysql_local_infile=1;' % (
                conf['db_server'],
                conf['db_port'],
                memory / 2,
                memory / 4,
                memory / 4)
        omcl_conf['dbLogin'] = conf['db_login']
        omcl_conf['dbPassword'] = conf['db_password']

        log.debug('Database connection string: ' + omcl_conf['dbConnectString'])

    with open(join(working_dir, orthomcl_config_fname), 'w') as ocf:
        ocf.writelines('='.join(item) + '\n' for item in omcl_conf.items())


def get_starting_step(start_from, log_file):
    if start_from == 'uselog':
        with log_file as f:
            for l in f.readlines():
                if 'Done' in l:
                    return None, l[41:].strip()
    try:
        return int(start_from), None
    except ValueError:
        return start_from, None