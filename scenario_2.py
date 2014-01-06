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
                 assembly, genes, proteome, existing_blast_result, existing_proteomes,
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

    result = workflow.run(start_after, start_from, overwrite, ask_before)
    if result == 0:
        log.info('Done.')
        log.info('Log in ' + join(working_dir, log_fname))
        log.info('Groups in ' + join(working_dir, steps.orthogroups_file))
    return result


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous genes.')

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

    indent = ' ' * len('usage: ' + basename(__file__) + ' ')
    op.usage = basename(__file__) + ' [--existing-blast-results TSV]\n' + \
        indent + '[--existing_proteomes DIR]\n' + \
        indent + '[--assembly FASTA]\n' + \
        indent + '[--genes GFF]\n' + \
        indent + '[--proteome FASTA]\n'

    add_common_arguments(op, indent)

    params = op.parse_args(args)

    if params.existing_blast_results and not params.existing_proteomes:
        interrupt('You need also provide existing proteomes.')

    if not params.species_file and \
       not params.ids_file and \
       not params.annotations_dir and \
       not params.start_from:
        interrupt('Either --species-file, --ids-file, --annotations-dir, '
                  'or --start-from has to be specified.')
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

    if p.existing_proteomes: p.existing_proteomes = join(getcwd(), p.existing_proteomes)
    if p.existing_blast_results: p.existing_blast_results = join(getcwd(), p.existing_blast_results)

    run_workflow(working_dir=p.out_dir,

                 assembly=p.plus_assembly, genes=p.plus_gff, proteome=p.plus_proteome,
                 existing_blast_result=p.existing_blast_results, existing_proteomes=p.existing_proteomes,

                 ask_before=p.ask_each_step,
                 start_after=start_after, start_from=start_from, overwrite=True,
                 threads=p.threads,
                 proxy=p.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])