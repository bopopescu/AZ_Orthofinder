from genericpath import isfile, isdir
from os.path import join
import sys
from src import config


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous.')
    op.add_argument('--species-file', dest='species_file',
                    help='File with a list of organism names as in Genbank. '
                         'For example, "Salmonella enterica subsp. enterica '
                         'serovar Typhi str. P-stx-12".'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('--ids-file', dest='ids_file',
                    help='File with reference ids to fetch from Genbank.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('--annotations-dir', dest='annotations_dir',
                    help='Directory with .gb files.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('--existing-blast-results', dest='existing_blast_results',
                    help='Existing proteomes directory.')

    op.add_argument('--existing-proteomes', dest='existing_proteomes',
                    help='Existing blast tsv.')

    op.add_argument('-a', '--assembly', dest='plus_assembly',
                    help='Fasta file with contigs to process with an existing blast tsv, '
                         'specified with --blast-results.')

    op.add_argument('-g', '--genes', dest='plus_gff',
                    help='Genes gff file to process with an existing blast tsv, '
                         'specified with --blast-results.')

    op.add_argument('-p', '--proteome', dest='plus_proteome',
                    help='Proteins fasta file to process with an existing blast tsv, '
                         'specified with --blast-results.')

    op.add_argument('-o', dest='out_dir', required=True,
                    help='The directory that will contain the resulting group.txt file, '
                         'as well as intermediate results.')

    # op.add_argument('-w', '--overwrite', dest='overwrite', action='store_true', default=False,
    #                 help='By default, the tool reuses existing intermediate results.'
    #                      'This option makes the tool overwrite any existing data.')

    op.add_argument('--start-from', dest='start_from', default=0,
                    help='Start from the specified step. '
                         'Either name or number (see log.txt) or "uselog".'
                         'If "uselog", the last "Done" record in log.txt will be searched.')

    op.add_argument('-t', '--threads', dest='threads', default=1,
                    help='Number of threads to run Blast.')

    op.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)

    op.add_argument('--ask', '--ask-each-step',
                    dest='ask_each_step', action='store_true', default=False,
                    help='Wait for user to press ke every time before proceed to next step.')

    op.add_argument('--proxy', dest='proxy', default=None, help='Proxy for FTP, for example: '
                                                                '--proxy 198.260.1.1:3333')
    usage = 'find_orthologs [--species-file FILE] \n' \
     '                      [--ids-file FILE] \n' \
     '                      [--annotations-dir DIR] \n' \
     ' \n' \
     '                      [--existing-blast-results TSV] \n' \
     '                      [--existing_proteomes DIR] \n' \
     '                      [--assembly FASTA] \n' \
     '                      [--genes GFF] \n' \
     '                      [--proteome FASTA] \n' \
     ' \n' \
     '                       -o OUTPUT_DIRECTORY \n' \
     '                      [--start-from STEP_NAME] \n' \
     '                      [-t threads_num] \n' \

    op.usage = usage

    params = op.parse_args(args)

    def interrupt(msg):
        print >> sys.stderr, msg
        exit(1)

    def check_file(fpath):
        if fpath and not isfile(fpath):
            interrupt('File ' + fpath + ' does not exist or is a directory.')

    def check_dir(dirpath):
        if dirpath and not isdir(dirpath):
            interrupt('Directory ' + dirpath + ' does not exist or is a file.')

    if params.existing_proteomes and not params.existing_blast_results:
        interrupt('You need also provide existing blast results tsv.')

    elif params.existing_blast_results and not params.existing_proteomes:
        interrupt('You need also provide existing proteomes.')

    elif params.existing_blast_results and params.existing_proteomes:
        if not params.plus_assembly and \
           not params.plus_gff and \
           not params.plus_proteome and \
           not params.start_from:
            interrupt('Either --assembly, --genes, --proteome, '
                      'or --start-from has to be specified.')
        check_file(params.existing_blast_results)
        check_file(params.existing_proteomes)
        check_file(params.plus_assembly)
        check_file(params.plus_gff)
        check_file(params.plus_proteome)

    else:
        if not params.species_file and \
           not params.ids_file and \
           not params.annotations_dir and \
           not params.start_from:
            interrupt('Either --species-file, --ids-file, --annotations-dir, '
                      'or --start-from has to be specified.')
        check_file(params.species_file)
        check_file(params.ids_file)
        check_dir(params.annotations_dir)

    if params.start_from == 'uselog':
        if not isfile(join(params.out_dir, config.log_fname)):
            print >> sys.stderr, 'No %s in %s. Either check your path, or ' \
                                 'change the --start-from option' % \
                                 (config.log_fname, params.out_dir)
            exit(1)

    return params