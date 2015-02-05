#!/usr/bin/env python
from shutil import move, copy
import sys
import logging
from os import chdir, mkdir, getcwd, listdir, makedirs
from os.path import join, exists, isdir, isfile, dirname, realpath, \
    basename, relpath, splitext, abspath, expanduser
from Bio import SeqIO

from src.fetch_annotations import fetch_annotations_for_species_from_ftp, \
    fetch_annotations_for_ids, fetch_annotations_species_name_entrez
from src.make_proteomes import adjust_proteomes, make_proteomes
from src.steps import check_results_existence
from src import steps
from src.argparse import ArgumentParser

from src.utils import make_workflow_id, read_list, set_up_config, \
    get_starting_step, interrupt, register_ctrl_c, \
    test_entrez_conn, test_blast_conn, test_ftp_conn, check_and_install_tools

from src.parse_args import arg_parse_error, check_file, check_dir, \
    add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow, Step

from src import config

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


#def run_workflow(working_dir,
#                 species_list, ids_list, annotations, proteomes, prot_id_field,
#                 min_length, max_percent_stop, evalue,
#                 ask_before=False,
#                 start_after=None, start_from=None, overwrite=True,
#                 threads=1,
#                 proxy=None,
#                 **kwargs):

def parse_args(args):
    op = ArgumentParser(description='Find groups of orthologous genes.')

    #op.add_argument(dest='directory')
    op.add_argument('-o', '--out', '--dir', dest='out', required=True)

    op.add_argument('-g', '--gbs', '--annotations', dest='annotations')
    op.add_argument('-p', '--proteins', '--proteomes', dest='proteomes')
    op.add_argument('--no-download', dest='download_anno', action='store_false', default=True)
    op.add_argument('-s', '--species', '--species-list', dest='species_list')
    op.add_argument('-i', '--ids', '--ids-list', dest='ids_list')

    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)

    op.usage = '''Finding orthogroups for a list of annotations / proteomes / ref ids / species.

Test runs:
    python scenario_1.py --ids test_input/ids.txt -o test_ids
    python scenario_1.py --proteomes test_input/proteins -o test_proteomes

Usage: %s [-p <proteomes dir>] [-a <.gb files dir>] [-i <gb ids file>] [-s <strain names file>]
                     [-o <dir>] [--jobs 30] [--start-from <step num>]
    -o  Output directory.
    -g  Directory with .gb files for references with annotations.
    -p  Directory with fasta (or faa, fa) files of protein sequences. If they
        are named by their reference ids (i.e. NC_005816.1.fasta), annotations
        will be downloaded from NCBI.
    -i  File with reference ids (will be fetched from NCBI).
    -s  File with a list of organism names as in Genbank.

    --prot-id-field
        When specifying proteomes, use this fasta id field number
        to retrieve protein ids (default if 1, like >NC_005816.1|NP_995567.1 ...).
    ''' % basename(__file__)

    #-a  --annotations-dir  Directory with .gb files.
    #-p  --proteomes-dir    Directory with fasta files of proteomes.
    #-i  --ids-list         File with reference ids (will be fetched from Genbank).
    #-s  --species-list     File with a list of organism names as in Genbank.
    #                       For example, "Salmonella enterica subsp. enterica serovar Typhi str. P-stx-12".
    #'''

    #indent = ' ' * len('usage: ' + basename(__file__) + ' ')
    #op.usage = basename(__file__) + ' [--annotations-dir DIR]\n' + \
    #    indent + '[--proteomes-dir DIR]\n' + \
    #    indent + '[--ids-file FILE]\n' + \
    #    indent + '[--species-file FILE]\n'

    add_common_arguments(op)

    p = op.parse_args(args)

    check_common_args(p)

    if not p.out:
        arg_parse_error('Specify output directory with -o.')
    if isfile(expanduser(p.out)):
        arg_parse_error('%s is a file' % p.out)
    p.out = abspath(expanduser(p.out))
    if not isdir(p.out):
        makedirs(p.out)

    if p.species_list:
        check_file(expanduser(p.species_list))
        p.species_list = abspath(expanduser(p.species_list))

    if p.ids_list:
        check_file(expanduser(p.ids_list))
        p.ids_list = abspath(expanduser(p.ids_list))

    if p.proteomes:
        check_dir(expanduser(p.proteomes))
        p.proteomes = abspath(expanduser(p.proteomes))

    if p.annotations:
        check_dir(expanduser(p.annotations))
        p.annotations = abspath(expanduser(p.annotations))

    #if p.species_list or p.ids_list:
    #    if not isdir(p.out):
    #        mkdir(p.out)
    #else:
    #    if not p.directory:
    #        arg_parse_error('Directory or file must be specified.')
    #        check_dir(p.directory)
    return p


def collect_proteomes_and_annotaitons(input_dir):
    proteomes = []
    annotations = []

    files = listdir(input_dir)
    if not files:
        interrupt('Directory contains no files.')

    for f in (join(input_dir, f) for f in files if isfile(join(input_dir, f))):
        if '.' in f and splitext(f)[1] in ['.fasta', '.faa', '.fa', '.fsa']:
            try:
                log.debug('   Checking if %s is fasta.' % f)
                next(SeqIO.parse(f, 'fasta'))
            except ValueError, e:
                pass
            else:
                proteomes.append(f)
                continue

        if '.' in f and splitext(f)[1] in ['.gb', '.genbank', '.gbk']:
            try:
                log.debug('   Checking if %s is genbank.' % f)
                SeqIO.read(f, 'genbank')
            except Exception, e:
                log.debug(str(e) + ', ' + f)
            else:
                annotations.append(f)

    log.debug('')
    return proteomes, annotations


def step_prepare_proteomes_and_annotations(p):
    def run(starting_from_here=False):
        if p.species_list:
            if not test_entrez_conn():
                log.error('   No internet connection: cannot fetch annotations.')
                return 4

            log.debug('   Using species list: ' + str(p.species_list))
            gb_ids = read_list(p.species_list)
            log.debug('species_list: ' + str(gb_ids))
            res = fetch_annotations_species_name_entrez(config.annotations_dir, gb_ids, p.proxy)
            if res != 0: return res
            return make_proteomes(config.annotations_dir, config.proteomes_dir)

        elif p.ids_list:
            if not test_entrez_conn():
                log.error('No internet connection: cannot fetch annotations.')
                return 4

            log.debug('   Using ref ids: ' + str(p.ids_list))
            ref_ids = read_list(p.ids_list)
            res = fetch_annotations_for_ids(config.annotations_dir, ref_ids)
            if res != 0: return res
            return make_proteomes(config.annotations_dir, config.proteomes_dir)

        else:
            proteomes, annotations = [], []

            if p.proteomes:
                proteomes, annotations = collect_proteomes_and_annotaitons(p.proteomes)
                if proteomes == []:
                    interrupt('No fasta found in ' + p.proteomes)

            if p.annotations:
                proteomes, annotations = collect_proteomes_and_annotaitons(p.annotations)
                if annotations == []:
                    interrupt('No gb files found in ' + p.annotations)

            #if not proteomes and not annotations:
            #    interrupt('Directory must contain fasta or genbank files.')
            #
            #if proteomes and annotations:
            #    log.warn('Directory %s contains both fasta and genbank files, using fasta.')

            if annotations:
                if not isdir(config.annotations_dir):
                    mkdir(config.annotations_dir)

                for annotation in annotations:
                    copy(annotation, config.annotations_dir)

                return make_proteomes(config.annotations_dir, config.proteomes_dir)

            elif proteomes:
                if not isdir(config.proteomes_dir):
                    mkdir(config.proteomes_dir)

                if p.download_anno:
                    if not test_entrez_conn():
                        #log.error('   Error: no internet connection, cannot fetch annotations. '
                        #          'You can start without a --no-fetch option, in this case '
                        #          'a reduced version of orthogroups.txt with no annotations will be produced.')
                        #return 1
                        log.error('   Warning: no internet connection, cannot fetch annotations. '
                                  'A reduced version of orthogroups.txt with no annotations will be produced.')
                    else:
                        # ref_ids = [splitext(basename(prot_file))[0] for prot_file in proteomes]
                        # fetch_annotations_for_ids(config.annotations_dir, ref_ids)

                        gb_ids = [splitext(basename(prot_file))[0] for prot_file in proteomes]
                        log.debug('ids_list: ' + str(gb_ids))
                        res = fetch_annotations_for_ids(config.annotations_dir, gb_ids, p.proxy)
                        if res > 0:
                            return res
                        if res == -1:
                            p.download_anno = False

                return adjust_proteomes(proteomes, config.proteomes_dir, p.prot_id_field)

    return Step(
        'Preparing proteomes and annotations',
        run=run,
        prod_files=[config.proteomes_dir, config.annotations_dir])


def main(args):
    register_ctrl_c()

    p = parse_args(args)
    log_path = set_up_logging(p.debug, p.out)
    log.info('python ' + __file__ + ' ' + ' '.join(args))
    log.info('logging to ' + log_path)
    log.info('')

    try:
        working_dir = p.out

        with open(config.config_file) as f:
            conf = dict(l.strip().lower().split('=', 1)
                        for l in f.readlines() if l.strip() and l.strip()[0] != '#')

        check_and_install_tools(p.debug, conf.get('db_vendor', 'sqlite') == 'sqlite', log_path)

        start_from, start_after = get_starting_step(p.start_from, join(p.out, log_fname))

        log.debug('Changing to %s' % working_dir)
        chdir(working_dir)

        if not p.overwrite:
            check_results_existence()

        if not exists(config.intermediate_dir):
            mkdir(config.intermediate_dir)

        set_up_config(working_dir)

        # Building the workflow
        workflow = Workflow(working_dir, id=make_workflow_id(working_dir),
                            cmdline_args=['python', __file__] + args)
        log.debug('Workflow id is "' + workflow.id + '"')
        log.debug('')

        suffix = '' if conf.get('db_vendor', 'sqlite') == 'sqlite' else '_' + workflow.id

        njobs = p.threads or p.jobs or 30

        workflow.extend([
            step_prepare_proteomes_and_annotations(p),
            steps.filter_proteomes(
                min_length=int(p.min_length),
                max_percent_stop=int(p.max_percent_stop)),
            steps.make_blast_db(),
            steps.blast(
                workflow.id,
                int(p.threads) or int(p.jobs) or 30,
                on_cluster=njobs and not p.threads,
                evalue=float(p.evalue)),
            steps.parse_blast_results(),
            steps.clean_database(suffix),
            steps.install_schema(suffix),
            steps.load_blast_results(suffix),
            steps.find_pairs(suffix),
            steps.dump_pairs_to_files(suffix),
            steps.mcl(p.debug),
            steps.step_save_orthogroups()])

        result = workflow.run(
            start_after, start_from,
            overwrite=True,
            ask_before=p.ask_each_step)

        if result == 0:
            log.info('Done.')
            log.info('Log is in ' + join(working_dir, log_fname))
            if isfile(join(working_dir, config.orthogroups_file)):
                log.info('Groups are in ' + join(working_dir, config.orthogroups_file))
                if isfile(config.nice_orthogroups_file):
                    log.info('Groups with aligned columns are in ' +
                             join(working_dir, config.nice_orthogroups_file))
            else:
                log.info('Groups in short format are in ' +
                         join(working_dir, config.short_orthogroups_file))

            if isfile(log_fname):
                with open(log_fname, 'a') as f:
                    f.write('\n')

        return result

    except (KeyboardInterrupt, SystemExit, GeneratorExit):
        if isfile(log_fname):
            with open(log_fname, 'a') as f:
                f.write('\n')
        return 1

    except Exception as e:
        log.error('')
        log.exception('Unexpected error!')
        if isfile(log_fname):
            with open(log_fname, 'a') as f:
                f.write('\n')
        return 2


if __name__ == '__main__':
    exit(main(sys.argv[1:]))