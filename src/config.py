from os.path import dirname, realpath, join, basename

BLAST_DBSIZE = 100000

config_file = join(dirname(realpath(__file__)), '../config.txt')

src_dir = dirname(realpath(__file__))
orthomcl_sqlite_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software/bin')
orthomcl_mysql_bin_dir = join(dirname(realpath(__file__)), 'orthomcl_software_mysql/bin')
orthomcl_bin_dir = orthomcl_sqlite_bin_dir
with open(config_file) as f:
    conf = dict(l.strip().lower().split('=', 1)
                for l in f.readlines() if l.strip() and l.strip()[0] != '#')
    if conf['db_vendor'] == 'sqlite':
        orthomcl_bin_dir = orthomcl_sqlite_bin_dir
    else:
        orthomcl_bin_dir = orthomcl_mysql_bin_dir
mcl_dir = join(dirname(realpath(__file__)), 'mcl_software')
mysql_linux_tar = join(dirname(realpath(__file__)), 'mysql-5.6.15-linux-glibc2.5-x86_64.tar.gz')
mysql_osx_tar = join(dirname(realpath(__file__)), 'mysql-5.6.15-osx10.7-x86_64.tar.gz')
mysql_extracted_dir = join(dirname(realpath(__file__)), 'mysql')
mysql_cnf = join(dirname(realpath(__file__)), 'mysql.cnf')
sqlite_fname = 'sqlite.db'

log_fname = 'log.txt'
logger_name = 'orthofinder'


proteomes_dir             = 'proteomes'
annotations_dir           = 'annotations'
intermediate_dir          = 'intermediate'
sqlite_file               = 'intermediate/sqlite.db'
sql_log                   = 'intermediate/log.sql'
good_proteins             = 'intermediate/good_proteins.fasta'
poor_proteins             = 'intermediate/poor_proteins.fasta'
blast_db                  = 'intermediate/blastdb'
blast_out                 = 'intermediate/blasted.tsv'
similar_sequences         = 'intermediate/similar_sequences.txt'
pairs_log                 = 'intermediate/orthomclpairs.log'
mcl_input                 = 'intermediate/mcl_input'
mcl_output                = 'intermediate/mcl_output'
pairs_dir                 = 'intermediate/pairs'
pairs_orthologs           = 'intermediate/pairs/potentialOrthologs.txt'
pairs_inparalogs          = 'intermediate/pairs/potentialInparalogs.txt'
pairs_coorthologs         = 'intermediate/pairs/potentialCoorthologs.txt'
groups_file               = 'groups.txt'
singletons_file           = 'singletons.txt'
orthogroups_file          = 'orthogroups.tsv'
nice_orthogroups_file     = 'orthogroups_nice.txt'
short_orthogroups_file    = 'orthogroups_short.tsv'
assembly_singletones_file = 'assembly_singletones.txt'
singletone_dir            = 'new_singletones'


orthomcl_config_fname = 'orthomcl.config'
orthomcl_config_final_path = join(intermediate_dir, orthomcl_config_fname)