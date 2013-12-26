from os.path import dirname, realpath, join, basename

BLAST_DBSIZE = 100000

orthomcl_config = join(dirname(realpath(__file__)), 'orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')
config_file = join(dirname(realpath(__file__)), '../config.txt')

log_fname = 'log.txt'
logger_name = 'orthofinder'
