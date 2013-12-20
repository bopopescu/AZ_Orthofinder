import logging
import sys
from os.path import dirname, realpath, join

BLAST_DBSIZE = 100000

orthomcl_config = join(dirname(realpath(__file__)), 'orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')
config = join(dirname(realpath(__file__)), '../config.txt')

log_fname = 'log.txt'
logger_name = 'orthofinder'

def set_up_logging(debug, working_dir):
    logger = logging.getLogger(log_fname)
    logger.setLevel(logging.DEBUG)

    class InfoFilter(logging.Filter):
        def filter(self, rec):
            return rec.levelno in (logging.DEBUG, logging.INFO)

    console_formatter = logging.Formatter(
        '%(asctime)-15s  %(message)s' if debug else '%(message)s',
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

    fh = logging.FileHandler(join(working_dir, log_fname), 'a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)-15s  %(levelname)-8s  %(message)s',
        datefmt='%c'))
    logger.addHandler(fh)