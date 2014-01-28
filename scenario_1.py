#!/usr/bin/env python
from shutil import move, copy
import sys
import logging
from os import chdir, mkdir, getcwd, listdir
from os.path import join, exists, isdir, isfile, dirname, realpath, \
    basename, relpath, splitext, abspath
from Bio import SeqIO

from src.fetch_annotations import fetch_annotations_for_species_from_ftp, \
    fetch_annotations_for_ids
from src.make_proteomes import adjust_proteomes, make_proteomes
from src.steps import check_results_existence
from src import steps
from src.argparse import ArgumentParser

from src.utils import make_workflow_id, read_list, set_up_config, \
    get_starting_step, interrupt, register_ctrl_c, test_internet_conn, \
    check_and_install_tools
from src.parse_args import arg_parse_error, check_file, check_dir, \
    add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow, Step

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

    op.add_argument('-g', '--annotations', '--gbs', dest='annotations')
    op.add_argument('-p', '--proteins', '--proteomes', dest='proteomes')
    op.add_argument('-s', '--species', '--species-list', dest='species_list')
    op.add_argument('-i', '--ids', '--ids-list', dest='ids_list')

    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)

    op.usage = '''Finding orthogroups for a list of annotations / proteomes / ref ids / species.

    usage: %s [--proteomes dir] [--annotations dir] [--ids-list file] [--species-list file]
                        [-o] [-t num] [--start-from step]

    -o:                  Output directory.

    Optional arguments:
    -g --gbs:            Directory with gb files.

    -p --proteomes:      Directory with fasta (or faa) protein files, named by their reference ids
                         (i.e. NC_005816.1.fasta). Can contain annotations from Prodigal.

    -s --species-list:   File with a list of organism names as in Genbank.

    -i --ids-list:       File with reference ids (will be fetched from Genbank).

    --prot-id-field:     When specifying proteomes, use this fasta id field number
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
    if isfile(p.out):
        arg_parse_error('%s is a file' % p.out)
    p.out = abspath(p.out)
    if not isdir(p.out):
        mkdir(p.out)

    if p.species_list:
        check_file(p.species_list)
        p.species_list = abspath(p.species_list)

    if p.ids_list:
        check_file(p.ids_list)
        p.ids_list = abspath(p.ids_list)

    if p.proteomes:
        check_dir(p.species_list)
        p.proteomes = abspath(p.proteomes)

    if p.annotations:
        check_dir(p.species_list)
        p.annotations = abspath(p.annotations)

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


def step_prepare_proteomes_and_annotations(p, internet_is_on):
    def run():
        if p.species_list:
            if not internet_is_on:
                log.error('   No internet connection: cannot fetch annotations.')
                return 4

            log.debug('   Using species list: ' + str(p.species_list))
            species_list = read_list(p.species_list)
            log.debug('species_list: ' + str(species_list))
            res = fetch_annotations_for_species_from_ftp(steps.annotations_dir, species_list, p.proxy)
            if res != 0: return res
            return make_proteomes(steps.annotations_dir, steps.proteomes_dir)

        elif p.ids_list:
            if not internet_is_on:
                log.error('No internet connection: cannot fetch annotations.')
                return 4

            log.debug('   Using ref ids: ' + str(p.ids_list))
            ref_ids = read_list(p.ids_list)
            res = fetch_annotations_for_ids(steps.annotations_dir, ref_ids)
            if res != 0: return res
            return make_proteomes(steps.annotations_dir, steps.proteomes_dir)

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
                if not isdir(steps.annotations_dir):
                    mkdir(steps.annotations_dir)

                for annotation in annotations:
                    copy(annotation, steps.annotations_dir)

                return make_proteomes(steps.annotations_dir, steps.proteomes_dir)

            elif proteomes:
                if not isdir(steps.proteomes_dir):
                    mkdir(steps.proteomes_dir)

                if not internet_is_on:
                    log.warn('   Warning: no internet connection, cannot fetch annotations. '
                             'A reduced version of orthogroups.txt with no annotations will be produced.')
                else:
                    ref_ids = [splitext(basename(prot_file))[0] for prot_file in proteomes]
                    fetch_annotations_for_ids(steps.annotations_dir, ref_ids)

                return adjust_proteomes(proteomes, steps.proteomes_dir, p.prot_id_field)

    return Step(
       'Preparing proteomes and annotations',
        run=run,
        prod_files=[steps.proteomes_dir, steps.annotations_dir])


def main(args):
    register_ctrl_c()

    p = parse_args(args)
    log_path = set_up_logging(p.debug, p.out)
    log.info('python ' + basename(__file__) + ' ' + ' '.join(args))
    log.info('logging to ' + log_path)
    log.info('')

    check_and_install_tools(p.debug, log_path)
    set_up_config()

    start_from, start_after = get_starting_step(p.start_from, join(p.out, log_fname))

    working_dir = p.out
    log.info('Changing to %s' % working_dir)
    chdir(working_dir)

    workflow = Workflow(working_dir, id=make_workflow_id(working_dir))
    log.info('Workflow id is "' + workflow.id + '"')
    log.info('')
    suffix = '_' + workflow.id

    if not p.overwrite:
        check_results_existence()

    if not exists(steps.intermediate_dir):
        mkdir(steps.intermediate_dir)

    internet_is_on = test_internet_conn()

    workflow.extend([
        step_prepare_proteomes_and_annotations(p, internet_is_on),

        steps.filter_proteomes(
            min_length=int(p.min_length),
            max_percent_stop=int(p.max_percent_stop)),
        steps.make_blast_db(),
        steps.blast(p.threads, evalue=float(p.evalue)),
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
        log.info('Groups are in ' + join(working_dir, steps.orthogroups_file))
        if isfile(steps.nice_orthogroups_file):
            log.info('Groups with aligned columns are in ' +
                     join(working_dir, steps.nice_orthogroups_file))
    return result


if __name__ == '__main__':
    exit(main(sys.argv[1:]))