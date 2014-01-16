from os import environ, access, X_OK, remove, getcwd, chdir, mkdir
from os.path import basename, isdir, split, join, pathsep, isfile, dirname, exists, devnull
from random import randint
from shutil import copy
from subprocess import call, Popen, PIPE
from config import config_file, orthomcl_config, log_fname
import logging

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

    return (basename(working_dir) or basename(dirname(working_dir))).\
        replace(' ', '_').replace('/', '_')


def interrupt(msg, code=1):
    log.error(msg)
    exit(code)


def register_ctrl_c():
    import signal

    def signal_handler(signal, frame):
        print ''
        exit(0)

    signal.signal(signal.SIGINT, signal_handler)


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


def check_installed_tools(tools):
    ok = True
    for tool in tools:
        if not which(tool):
            log.warn('"' + tool + '" might not be installed.')
            ok = False
    if not ok:
        exit(3)


def check_and_install_mcl(mcl_src_path, log_file_path):
    if which('mcl'):
        return 'mcl'

    mcl_bin_path = join(mcl_src_path, 'bin', 'bin', 'mcl')
    if not exists(mcl_bin_path):
        with open(log_file_path) as h:
            log.info('Compiling MCL...')
            cur_dir = getcwd()
            chdir(mcl_src_path)
            mcl_path = join(mcl_src_path, 'bin')

            def run(command, exit_code=1):
                log.info('   ' + command)
                if call(command.split(), stderr=h, stdout=open(devnull, 'w')) != 0:
                    log.error('Cannot install mcl :( Try manually.')
                    exit(2)

            run('./configure -q --prefix=' + mcl_path, 1)
            run('make', 2)
            run('make check', 3)
            run('make install', 4)

        chdir(cur_dir)
    return mcl_bin_path


def check_perl_modules(src_path, log_file_path, debug):
    def check(env):
        dbimysql = call('perl -MDBD::mysql -e 1'.split(),
                         stderr=open(devnull, 'w'),
                         stdout=open(devnull, 'w'),
                         env=env) == 0
        if not dbimysql:
            log.debug('DBD::mysql is not installed.')

        dbd = call('perl -MDBI -e 1'.split(),
                    stderr=open(devnull, 'w'),
                    stdout=open(devnull, 'w'),
                    env=env) == 0
        if not dbd:
            log.debug('DBI is not installed.')

        return dbimysql and dbd

    perl_modules_path = join(src_path, 'perl_modules')
    lib_path = join(perl_modules_path, 'lib')
    man_path = join(perl_modules_path, 'man')
    cpan_path = join(src_path, '.cpan')
    cpan_cpan_path = join(cpan_path, 'CPAN')

    if not isdir(perl_modules_path):
        mkdir(perl_modules_path)
        mkdir(lib_path)
        mkdir(man_path)
        mkdir(join(man_path, 'man1'))
        mkdir(join(man_path, 'man3'))
    if not isdir(cpan_path): mkdir(cpan_path)
    if not isdir(cpan_cpan_path): mkdir(cpan_cpan_path)

    with open(join(cpan_cpan_path, 'MyConfig.pm'), 'w') as cpan_f, \
         open(join(cpan_cpan_path, 'MyConfig_template.pm')) as cpan_tmpl_f:
        cpan_f.write(cpan_tmpl_f.read().replace('{PATH}', src_path))

    if 'PERL5LIB' not in environ:
        environ['PERL5LIB'] = lib_path
    else:
        environ['PERL5LIB'] = environ['PERL5LIB'] + ':' + lib_path

    if 'MANPATH' not in environ:
        environ['MANPATH'] = lib_path
    else:
        environ['MANPATH'] = environ['MANPATH'] + ':' + man_path

    env = {'PERL5LIB': lib_path, 'MANPATH': man_path}
    if check(environ):
        return

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

    log.info('Installing DBD::mysql...')
    cmd = ['cpan', '-f', '-j', join(cpan_cpan_path, 'MyConfig.pm'), 'DBD::mysql']
    log.debug(' '.join(cmd))
    call(cmd, env=env)


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

    if not check(env):
        log.error('Could not install Perl modules. ')
        print '''Pleasy, try the following manually:
    $ perl -MCPAN -e shell
    cpan> o conf makepl_arg "mysql_config=mysql_config"
    cpan> install Data::Dumper
    cpan> install DBI
    cpan> force install DBD::mysql
    '''
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


def test_internet_conn():
    import urllib2
    try:
        response = urllib2.urlopen('http://74.125.228.100', timeout=1)
    except urllib2.URLError as err:
        return False
    else:
        return True


def set_up_config():
    with open(config_file) as cf:
        conf = dict(l.strip().split('=', 1) for l
                    in cf.readlines() if l.strip()[0] != '#')
        log.debug('Read conf: ' + str(conf))

    with open(orthomcl_config) as ocf:
        omcl_conf = dict(l.strip().split('=', 1) for l
                         in ocf.readlines() if l.strip()[0] != '#')

    omcl_conf['dbConnectString'] = \
        'dbi:mysql:database=orthomcl;host=127.0.0.1;port=%s;mysql_local_infile=1' % \
            conf['db_port']
    omcl_conf['dbLogin'] = conf['db_login']
    omcl_conf['dbPassword'] = conf['db_password']

    with open(orthomcl_config, 'w') as ocf:
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