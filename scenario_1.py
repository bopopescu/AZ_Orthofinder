#!/usr/bin/env python
from genericpath import isfile

import sys
import logging
from os import chdir, mkdir, getcwd, listdir
from os.path import join, exists, isdir, dirname, realpath, basename, relpath, splitext
from Bio import SeqIO
from src import steps
from src.argparse import ArgumentParser

from src.utils import which, make_workflow_id, read_list, set_up_config, get_start_after_from
from src.parse_args import interrupt, check_file, check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


def run_workflow(working_dir,
                 species_list, ids_list, annotations, proteomes, prot_id_field,
                 min_length, max_percent_stop,
                 ask_before=False,
                 start_after=None, start_from=None, overwrite=True,
                 threads=1,
                 proxy=None,
                 **kwargs):

    if not which('blastp') or not which('mcl'):
        if not which('blastp'):
            log.error('blastp installation required.')
        if not which('mcl'):
            log.error('mcl installation required.')
        return 3

    log.info('Changing to %s' % working_dir)
    chdir(working_dir)

    if not exists('intermediate'): mkdir('intermediate')

    workflow = Workflow(working_dir, id=make_workflow_id(working_dir))
    log.info('Workflow id is "' + workflow.id + '"')
    log.info('')
    suffix = '_' + workflow.id

    if species_list:
        log.debug('Using species list: ' + str(species_list))
        workflow.add(steps.step_fetching_annotations_for_species(species_list, proxy))
        workflow.add(steps.step_make_proteomes())

    elif ids_list:
        log.debug('Using ref ids: ' + str(ids_list))
        workflow.add(steps.step_fetch_annotations_for_ids(ids_list))
        workflow.add(steps.step_make_proteomes())

    elif annotations:
        log.debug('Using genbank files.')
        workflow.add(steps.step_make_proteomes(annotations))

    elif proteomes:
        log.debug('Using proteomes_files.')
        workflow.add(steps.step_adjust_proteomes(proteomes, prot_id_field))

    elif start_from == 0:
        log.error('Either species names, reference ids, annotations, proteomes, '
                  'or step to start from has to be be speciefied.')
        exit(1)

    workflow.extend([
        steps.filter_proteomes(min_length, max_percent_stop),
        steps.make_blast_db(),
        steps.blast(threads),
        steps.parse_blast_results(),
        steps.clean_database(suffix),
        steps.install_schema(suffix),
        steps.load_blast_results(suffix),
        steps.find_pairs(suffix),
        steps.dump_pairs_to_files(suffix),
        steps.mcl(),
        steps.step_save_orthogroups(annotations or None)])

    result = workflow.run(start_after, start_from, overwrite, ask_before)
    if result == 0:
        log.info('Done.')
        log.info('Log is in ' + join(working_dir, log_fname))
        log.info('Groups are in ' + join(working_dir, steps.orthogroups_file))
        log.info('Groups with aligned columns are in ' + join(working_dir, steps.nice_orthogroups_file))
    return result


def parse_args(args):
    op = ArgumentParser(description='Find groups of orthologous genes.')

    op.add_argument(dest='directory')
    op.add_argument('-s', '--species-list', dest='species_list')
    op.add_argument('-i', '--ids-list', dest='ids_list')
    op.add_argument('--annotations-files', dest='annotations_files')
    op.add_argument('--proteomes-files', dest='proteomes_files')
    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)

    op.usage = '''Finding orthogroups for a list of annotations / proteomes / ref ids / species.

    usage: %s directory [-t num] [--start-from step] [-i file] [-s file]

    Directory contains fasta files with proteomes.
    The directory can alternatively contain .gb files, or you can pass
    a file instead with list of reference ids or species: annotations
    will be fetched from Genbank instead.

    Optional arguments:
    -s:                File with a list of organism names as in Genbank.

    -i:                File with reference ids (will be fetched from Genbank).

    --prot-id-field:   When specifying proteomes, use this fasta id field number
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

    params = op.parse_args(args)

    check_common_args(params)

    if params.species_list or params.ids_list:
        if not isdir(params.directory): mkdir(params.directory)
        if params.species_list: check_file(params.species_list)
        if params.ids_list: check_file(params.ids_list)
        return params

    else:
        if not params.directory:
            interrupt('Directory or file must be specified.')
        return params


def main(args):
    p = parse_args(args)

    annotations = []
    proteomes = []
    species_list = []
    ref_id_list = []

    if not isdir(p.directory):
        interrupt('No such directory: ' + p.directory)

    set_up_logging(p.debug, p.directory)
    log.info(basename(__file__) + ' ' + ' '.join(args) + '\n')
    set_up_config()
    start_from, start_after = get_start_after_from(p.start_from, join(p.directory, log_fname))

    if p.species_list or p.ids_list:
        species_list = read_list(p.species_list, p.directory)
        ref_id_list = read_list(p.ids_list, p.directory)

    else:
        if 'proteomes' in listdir(p.directory):
            proteomes = [relpath(join('proteomes', f))
                         for f in listdir(join(p.directory, 'proteomes'))
                         if f and f[0] != '.']

        if 'annotations' in listdir(p.directory):
            annotations = [relpath(join('annotations', f))
                           for f in listdir(join(p.directory, 'annotations'))
                           if f and f[0] != '.']

        if not proteomes and not annotations:
            files = listdir(p.directory)
            if not files: interrupt('Directory contains no files.')

            for f in (join(p.directory, f) for f in files if isfile(join(p.directory, f))):
                if '.' in f and splitext(f)[1] in ['.fasta', '.faa', '.fa', '.fsa']:
                    try:
                        log.debug('   Checking if %s is fasta.' % f)
                        next(SeqIO.parse(f, 'fasta'))
                    except Exception, e:
                        pass
                    else:
                        proteomes.append(relpath(f, p.directory))
                        continue
                if '.' in f and splitext(f)[1] in ['.gb', '.genbank', '.gbk']:
                    try:
                        log.debug('   Checking if %s is genbank.' % f)
                        SeqIO.read(f, 'genbank')
                    except Exception, e:
                        log.debug(str(e) + ', ' + f)
                    else:
                        annotations.append(relpath(f, p.directory))
            log.debug('')

        if not proteomes and not annotations:
            interrupt('Directory must contain fasta or genbank files.')

        if proteomes and annotations:
            log.warn('Directory %s contains both fasta and genbank files, using fasta.')

        #if annotations: annotations = [join(getcwd(), path) for path in annotations]
        #if proteomes: proteomes = [join(getcwd(), path) for path in proteomes]

    return run_workflow(
        working_dir=p.directory,

        species_list=species_list, ids_list=ref_id_list,
        annotations=annotations, proteomes=proteomes,
        prot_id_field=int(p.prot_id_field),
        min_length=int(p.min_length), max_percent_stop=int(p.max_percent_stop),

        ask_before=p.ask_each_step,
        start_after=start_after, start_from=start_from, overwrite=True,
        threads=p.threads,
        proxy=p.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])