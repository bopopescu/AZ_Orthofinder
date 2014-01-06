from genericpath import isfile, isdir
from os.path import join
import sys
from src import config


def add_common_arguments(op, indent):
    op.add_argument('-o', dest='out_dir', required=True,
                    help='The directory that will contain the resulting group.txt file, '
                         'as well as intermediate results.')

    op.add_argument('--start-from', dest='start_from', default=0,
                    help='Start from the specified step. '
                         'Either name or number (see log.txt) or "uselog".'
                         'If "uselog", the last "Done" record in log.txt will be searched.')

    op.add_argument('--ask', '--ask-each-step',
                    dest='ask_each_step', action='store_true', default=False,
                    help='Wait for user to press ke every time before proceed to next step.')

    op.add_argument('-t', '--threads', dest='threads', default=30,
                    help='Number of threads to run Blast.')

    op.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)

    # op.add_argument('-w', '--overwrite', dest='overwrite', action='store_true', default=False,
    #                 help='By default, the tool reuses existing intermediate results.'
    #                      'This option makes the tool overwrite any existing data.')

    op.add_argument('--proxy', dest='proxy', default=None, help='Proxy for FTP, for example: '
                                                                '--proxy 198.260.1.1:3333')
    op.usage += \
        indent + ' -o DIR\n' + \
        indent + '[--start-from STEP_NAME]\n' + \
        indent + '[-t threads_num]\n'


def check_common_args(params):
    if params.start_from == 'uselog':
        if not isfile(join(params.out_dir, config.log_fname)):
            interrupt('No %s in %s. Either check your path, or '
                      'change the --start-from option' %
                      (config.log_fname, params.out_dir))

def interrupt(msg):
    print >> sys.stderr, msg
    exit(1)

def check_file(fpath):
    if fpath and not isfile(fpath):
        interrupt('File ' + fpath + ' does not exist or is a directory.')

def check_dir(dirpath):
    if dirpath and not isdir(dirpath):
        interrupt('Directory ' + dirpath + ' does not exist or is a file.')

