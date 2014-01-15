#!/usr/bin/env python
from genericpath import isfile

import sys
import logging
from os import chdir, mkdir, getcwd
from os.path import join, exists, isdir, dirname, realpath, basename, splitext, abspath
from src.make_proteomes import adjust_proteomes
from src.steps import check_results_existence
from src import steps

from src.utils import which, make_workflow_id, read_list, set_up_config, get_starting_step, register_ctrl_c, \
    check_installed_tools, test_internet_conn
from src.parse_args import arg_parse_error, check_file, check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow, Step, cmdline

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


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
    p = op.parse_args(args)

    check_common_args(p)

    if p.assembly: p.assembly = abspath(p.assembly)

    if not isdir(p.directory):
        arg_parse_error('Directory %s does not exist.' % p.directory)

    return p


def step_filter_assembly(assembly, assembly_name):
    return Step(
        'Filtering assembly',
         run=lambda: steps.filter_assembly(assembly, join(steps.intermediate_dir, assembly_name + '.fna')),
         req_files=[assembly],
         prod_files=[join(steps.intermediate_dir, assembly_name + '.fna')])


def finding_genes(assembly_name):
    return Step(
        'Finding genes',
         run=cmdline('prodigal',
             parameters=[
                 '-a', join(steps.proteomes_dir, assembly_name + '.fasta'),
                 '-o', join(steps.intermediate_dir, assembly_name + '.gff'),
                 '-i', join(steps.intermediate_dir, assembly_name + '.fna')]),
         req_files=[steps.proteomes_dir,
                    join(steps.intermediate_dir, assembly_name + '.fna')],
         prod_files=[join(steps.proteomes_dir, assembly_name + '.fasta')])


def step_adjust_new_proteome(assembly_name, id_field=0):
    return Step(
        'Adjusting proteome',
         run=lambda: adjust_proteomes(
             [join(steps.proteomes_dir, assembly_name + '.fasta')],
             steps.proteomes_dir,
             id_field),
         req_files=[steps.proteomes_dir,
                    join(steps.proteomes_dir, assembly_name + '.fasta')])


def step_blast_singletones():
    def blast_singletones(singletones_file, assembly_proteome):
        with open(singletones_file) as f:
            for line in f:
                strain, prot_id, gene_id, locus, product, description = line.split()
                # TODO: Blast!


    return Step(
        'Blasting singletones',
         run=lambda: blast_singletones(steps.assembly_singletones),
         req_files=[steps.assembly_singletones])


def main(args):
    register_ctrl_c()

    p = parse_args(args)
    set_up_logging(p.debug, p.directory)
    log.info('python ' + basename(__file__) + ' ' + ' '.join(args) + '\n')
    check_installed_tools(['blastp', 'mcl'])
    set_up_config()
    start_from, start_after = get_starting_step(p.start_from, join(p.directory, log_fname))

    working_dir = p.directory
    log.info('Changing to %s' % working_dir)
    chdir(working_dir)

    if not exists('intermediate'):
        arg_parse_error('You need to run scenario_1 on this directory first.')

    workflow = Workflow(working_dir, id=make_workflow_id(working_dir))
    log.info('Workflow id is "' + workflow.id + '"')
    log.info('')
    suffix = '_' + workflow.id

    internet_is_on = test_internet_conn()

    assembly_name = splitext(basename(p.assembly))[0]

    workflow.extend([
        step_filter_assembly(p.assembly, assembly_name),
        finding_genes(assembly_name),
        step_adjust_new_proteome(assembly_name, id_field=0),

        steps.filter_proteomes(min_length=int(p.min_length),
                               max_percent_stop=int(p.max_percent_stop)),
        steps.make_blast_db(),
        steps.blast(p.threads, assembly_name, evalue=float(p.evalue)),
        steps.parse_blast_results(),
        steps.clean_database(suffix),
        steps.install_schema(suffix),
        steps.load_blast_results(suffix),
        steps.find_pairs(suffix),
        steps.dump_pairs_to_files(suffix),
        steps.mcl(),
        steps.step_save_orthogroups(join(steps.proteomes_dir, assembly_name + '.fasta')),
        step_blast_singletones(),
    ])

    result = workflow.run(start_after, start_from, overwrite=True, ask_before=p.ask_each_step)
    if result == 0:
        log.info('Done.')
        log.info('Log is in ' + join(working_dir, log_fname))
        log.info('Groups are in ' + join(working_dir, steps.orthogroups_file))
        if isfile(steps.nice_orthogroups_file):
            log.info('Groups with aligned columns are in ' + join(working_dir, steps.nice_orthogroups_file))
    return result



if __name__ == '__main__':
    main(sys.argv[1:])