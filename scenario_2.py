#!/usr/bin/env python
from genericpath import isfile

import sys
import logging
from os import chdir, mkdir, getcwd
from os.path import join, exists, isdir, dirname, realpath, basename, splitext
from src import steps

from src.utils import which, make_workflow_id, read_list, set_up_config, get_start_after_from
from src.parse_args import interrupt, check_file, check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


def run_workflow(working_dir, assembly,
                 min_length, max_percent_stop, evalue,
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

    if not exists('intermediate'):
        interrupt('You need to run scenario_1 on this directory first.')

    workflow = Workflow(working_dir, id=make_workflow_id(working_dir))
    log.info('Workflow id is "' + workflow.id + '"')
    log.info('')
    suffix = '_' + workflow.id

    #if existing_blast_result and existing_proteomes:
    #    if assembly or genes:
    #        workflow.add(steps.step_predict_genes)
    #    if genes:
    #        workflow.add(steps.prepare_proteomes)
    #    steps.append(steps.step_blast_extend)

    assembly_name = splitext(basename(assembly))[0]

    workflow.extend([
        steps.step_filter_assembly(assembly, assembly_name),
        steps.finding_genes(assembly_name),
        steps.step_adjust_new_proteome(assembly_name, id_field=0),
        steps.filter_proteomes(),
        steps.make_blast_db(),
        steps.blast(threads, assembly_name),
        steps.parse_blast_results(),
        steps.clean_database(suffix),
        steps.install_schema(suffix),
        steps.load_blast_results(suffix),
        steps.find_pairs(suffix),
        steps.dump_pairs_to_files(suffix),
        steps.mcl(),
        steps.groups_to_files(workflow.id),
        steps.step_save_orthogroups()])

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

    op.add_argument(dest='directory')
    op.add_argument(dest='assembly')

    op.usage = '''Adding assembly to orthogroups.

    usage: %s directory assembly.fasta
    ''' % basename(__file__)

    #indent = ' ' * len('usage: ' + basename(__file__) + ' ')
    #op.usage = basename(__file__) + ' [--existing-blast-results TSV]\n' + \
    #    indent + '[--existing_proteomes DIR]\n' + \
    #    indent + '[--assembly FASTA]\n' + \
    #    indent + '[--genes GFF]\n' + \
    #    indent + '[--proteome FASTA]\n'

    add_common_arguments(op)

    params = op.parse_args(args)

    if not isdir(params.directory):
        interrupt('Directory %s does not exist.' % params.directory)

    check_common_args(params)

    return params


def main(args):
    p = parse_args(args)

    set_up_logging(p.debug, p.directory)
    log.info(basename(__file__) + ' ' + ' '.join(args) + '\n')
    set_up_config()
    start_from, start_after = get_start_after_from(p.start_from, join(p.directory, log_fname))

    p.assembly = join(getcwd(), p.assembly)

    run_workflow(working_dir=p.directory,

                 assembly=p.assembly,
                 #genes=p.plus_gff,
                 #proteome=p.plus_proteome,

                 min_length=int(p.min_length),
                 max_percent_stop=int(p.max_percent_stop),
                 evalue=float(p.evalue),

                 ask_before=p.ask_each_step,
                 start_after=start_after, start_from=start_from, overwrite=True,
                 threads=p.threads,
                 proxy=p.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])