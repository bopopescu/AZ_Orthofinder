import logging
import sys
from os.path import join

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

    std = logging.StreamHandler(sys.stdout)
    std.setLevel(logging.DEBUG if debug else logging.INFO)
    std.addFilter(InfoFilter())
    std.setFormatter(console_formatter)
    logger.addHandler(std)

    err = logging.StreamHandler(sys.stderr)
    err.setLevel(logging.WARN)
    err.setFormatter(console_formatter)
    logger.addHandler(err)

    log_fpath = join(working_dir, log_fname)
    fh = logging.FileHandler(log_fpath, mode)
    fh.setLevel(logging.DEBUG if debug else logging.INFO)
    fh.setFormatter(logging.Formatter(
        '%(asctime)-15s  %(message)s',
        datefmt='%c'))
    logger.addHandler(fh)

    with open(log_fpath, 'a') as f:
        f.write('*' * 50)

    return log_fpath
