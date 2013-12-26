#!/usr/bin/env python

import sys
import logging
from os import chdir, mkdir, getcwd
from os.path import join, exists, isdir, dirname, realpath
from src import steps

from src.utils import which, make_workflow_id, read_list, set_up_config, get_start_after_from
from src.parse_args import parse_args
from src.logger import set_up_logging
from src.Workflow import Workflow

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


def workflow_ini(working_dir,
                 assembly, genes, proteome, existing_blast_result, existing_proteomes,
                 species_list, ids_list, user_annotations_dir,
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

    if existing_blast_result and existing_proteomes:
        if assembly or genes:
            workflow.add(steps.step_predict_genes)
        if genes:
            workflow.add(steps.prepare_proteomes)
        steps.append(steps.step_blast_extend)

    else:
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

    if user_annotations_dir:
        workflow.add(steps.step_save_orthogroups(user_annotations_dir))
    else:
        workflow.add(steps.step_save_orthogroups())

    result = workflow.run(start_after, start_from, overwrite, ask_before)
    if result == 0:
        log.info('Done.')
        log.info('Log in ' + join(working_dir, log_fname))
        log.info('Groups in ' + join(working_dir, steps.orthogroups_file))
    return result


def main(args):
    p = parse_args(args[1:])

    if not isdir(p.out_dir): mkdir(p.out_dir)

    set_up_logging(p.debug, p.out_dir)

    log.info(' '.join(args) + '\n')

    set_up_config()

    start_from, start_after = get_start_after_from(p.start_from, join(p.out_dir, log_fname))

    species_list = read_list(p.species_file, p.out_dir)
    ref_id_list = read_list(p.ids_file, p.out_dir)
    if p.annotations_dir: p.annotations_dir = join(getcwd(), p.annotations_dir)
    if p.existing_proteomes: p.existing_proteomes = join(getcwd(), p.existing_proteomes)
    if p.existing_blast_results: p.existing_blast_results = join(getcwd(), p.existing_blast_results)

    workflow_ini(working_dir=p.out_dir,
                 assembly=p.plus_assembly, genes=p.plus_gff, proteome=p.plus_proteome,
                 existing_blast_result=p.existing_blast_results, existing_proteomes=p.existing_proteomes,
                 species_list=species_list, ids_list=ref_id_list, user_annotations_dir=p.annotations_dir,
                 ask_before=p.ask_each_step,
                 start_after=start_after, start_from=start_from, overwrite=True,
                 threads=p.threads,
                 proxy=p.proxy)


if __name__ == '__main__':
    main(sys.argv)