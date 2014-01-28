from os.path import dirname, realpath, join, basename

BLAST_DBSIZE = 100000

config_file = join(dirname(realpath(__file__)), '../config.txt')

src_dir = dirname(realpath(__file__))
orthomcl_config = join(dirname(realpath(__file__)), 'orthomcl.config')
orthomcl_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')
mcl_dir = join(dirname(realpath(__file__)), 'mcl_software')
mysql_linux_tar = join(dirname(realpath(__file__)), 'mysql-5.6.15-linux-glibc2.5-x86_64.tar.gz')
mysql_osx_tar = join(dirname(realpath(__file__)), 'mysql-5.6.15-osx10.7-x86_64.tar.gz')
mysql_extracted_dir = join(dirname(realpath(__file__)), 'mysql')
mysql_cnf = join(dirname(realpath(__file__)), 'mysql.cnf')

log_fname = 'log.txt'
logger_name = 'orthofinder'
