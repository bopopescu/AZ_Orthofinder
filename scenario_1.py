#!/usr/bin/env python
from genericpath import isfile

import sys
import logging
from os import chdir, mkdir, getcwd
from os.path import join, exists, isdir, dirname, realpath, basename
from src import steps

from src.utils import which, make_workflow_id, read_list, set_up_config, get_start_after_from
from src.parse_args import interrupt, check_file, check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


def run_workflow(working_dir,
                 species_list, ids_list, user_annotations_dir, proteomes_dir,
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
    log.info('Workflow id is "' + workflow.id + '"\n')
    suffix = '_' + workflow.id[:10] if workflow.id else ''

    if species_list:
        log.debug('Using species list: ' + str(species_list))
        workflow.add(steps.step_fetching_annotations_for_species(species_list, proxy))
        workflow.add(steps.step_make_proteomes())

    elif ids_list:
        log.debug('Using ref ids: ' + str(ids_list))
        workflow.add(steps.step_fetch_annotations_for_ids(ids_list))
        workflow.add(steps.step_make_proteomes())

    elif user_annotations_dir:
        log.debug('Using user_annotations_dir: ' + user_annotations_dir)
        workflow.add(steps.step_make_proteomes(user_annotations_dir))

    elif start_from == 0:
        log.error('Either species names, reference ids, annotations directory, '
                  'or step to start from has to be be speciefied.')
        exit(1)

    workflow.extend([
        steps.filter_proteomes(),
        steps.make_blast_db(),
        steps.blast(threads),
        steps.parse_blast_restults(),
        steps.clean_database(suffix),
        steps.install_schema(suffix),
        steps.load_blast_results(suffix),
        steps.find_pairs(suffix),
        steps.dump_pairs_to_files(suffix),
        steps.mcl()])

    workflow.add(steps.step_save_orthogroups(user_annotations_dir or None))

    result = workflow.run(start_after, start_from, overwrite, ask_before)
    if result == 0:
        log.info('Done.')
        log.info('Log is in ' + join(working_dir, log_fname))
        log.info('Groups are in ' + join(working_dir, steps.orthogroups_file))
        log.info('Groups with aligned columns are in ' + join(working_dir, steps.nice_orthogroups_file))
    return result


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous genes.')

    op.add_argument('-p', '--proteomes-dir', dest='proteomes_dir',
                    help='Directory with proteomes.')

    op.add_argument('-s', '--species-file', dest='species_file',
                    help='File with a list of organism names as in Genbank. '
                         'For example, "Salmonella enterica subsp. enterica '
                         'serovar Typhi str. P-stx-12".'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('-i', '--ids-file', dest='ids_file',
                    help='File with reference ids to fetch from Genbank.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    op.add_argument('-a', '--annotations-dir', dest='annotations_dir',
                    help='Directory with .gb files.'
                         'Either species-file, ids-file or annotations-dir'
                         'must be specified.')

    indent = ' ' * len('usage: ' + basename(__file__) + ' ')
    op.usage = basename(__file__) + ' [--proteomes-dir DIR]\n' + \
        indent + '[--ids-file FILE]\n' + \
        indent + '[--annotations-dir DIR]\n' + \
        indent + '[--species-file FILE]\n'

    add_common_arguments(op, indent)

    params = op.parse_args(args)

    if not params.proteomes_dir and \
       not params.species_file and \
       not params.ids_file and \
       not params.annotations_dir and \
       not params.start_from:
        interrupt('Either --proteomes-dir --species-file, --ids-file, --annotations-dir, '
                  'or --start-from has to be specified.')
    check_dir(params.proteomes_dir)
    check_file(params.species_file)
    check_file(params.ids_file)
    check_dir(params.annotations_dir)

    check_common_args(params)

    return params


def main(args):
    p = parse_args(args)
    if not isdir(p.out_dir): mkdir(p.out_dir)
    set_up_logging(p.debug, p.out_dir)
    log.info(' '.join(args) + '\n')
    set_up_config()
    start_from, start_after = get_start_after_from(p.start_from, join(p.out_dir, log_fname))

    species_list = read_list(p.species_file, p.out_dir)
    ref_id_list = read_list(p.ids_file, p.out_dir)
    if p.annotations_dir: p.annotations_dir = join(getcwd(), p.annotations_dir)
    if p.proteomes_dir: p.proteomes_dir = join(getcwd(), p.proteomes_dir)

    run_workflow(working_dir=p.out_dir,

                 species_list=species_list, ids_list=ref_id_list,
                 user_annotations_dir=p.annotations_dir, proteomes_dir=p.proteomes_dir,

                 ask_before=p.ask_each_step,
                 start_after=start_after, start_from=start_from, overwrite=True,
                 threads=p.threads,
                 proxy=p.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])