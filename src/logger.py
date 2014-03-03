import logging
import sys
from os.path import join
import time
import datetime

from config import log_fname


def set_up_logging(debug, working_dir, mode='a'):
    logger = logging.getLogger(log_fname)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    class InfoFilter(logging.Filter):
        def filter(self, rec):
            return rec.levelno in (logging.DEBUG, logging.INFO)

    console_formatter = logging.Formatter(
        '%(asctime)-15s  %(message)s',
        datefmt='%c')

    log_fpath = join(working_dir, log_fname)
    fh = logging.FileHandler(log_fpath, mode)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)-15s  %(message)s',
        datefmt='%c'))
    logger.addHandler(fh)

    std = logging.StreamHandler(sys.stdout)
    std.setLevel(logging.DEBUG if debug else logging.INFO)
    std.addFilter(InfoFilter())
    std.setFormatter(console_formatter)
    logger.addHandler(std)

    err = logging.StreamHandler(sys.stderr)
    err.setLevel(logging.WARN)
    err.setFormatter(console_formatter)
    logger.addHandler(err)

    with open(join(working_dir, log_fpath), 'a') as f:
        f.write('\n')
        f.write('*' * 30)
        today = datetime.datetime.now()
        f.write(' ' + today.strftime('%c') + ' ')
        f.write('*' * 30)
        f.write('\n\n')

    return log_fpath


def add_file_handler(working_dir, mode='a'):
    logger = logging.getLogger(log_fname)
    log_fpath = join(working_dir, log_fname)
    fh = logging.FileHandler(log_fpath, mode)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)-15s  %(message)s',
        datefmt='%c'))
    logger.addHandler(fh)
