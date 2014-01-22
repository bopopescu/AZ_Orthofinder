#!/usr/bin/env python
from genericpath import isfile
from shutil import copyfile, rmtree, copy

import sys
import logging
from os import chdir, mkdir, getcwd, listdir, symlink
from os.path import join, exists, isdir, dirname, realpath,\
    basename, splitext, abspath
from src.process_assembly import filter_assembly
from src.make_proteomes import adjust_proteomes
from src.steps import check_results_existence
from src import steps

from src.utils import which, make_workflow_id, read_list, \
    set_up_config, get_starting_step, register_ctrl_c, \
    check_installed_tools, test_internet_conn, \
    check_install_mcl, check_perl_modules, check_and_install_tools
from src.parse_args import arg_parse_error, check_file,\
    check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging
from src.Workflow import Workflow, Step, cmdline

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous genes.')

    op.add_argument(dest='directory')

    op.add_argument('-a', '--assemblies', dest='assemblies')
    op.add_argument('-g', '--gbs', dest='annotations')
    op.add_argument('-p', '--proteomes', dest='proteomes')
    op.add_argument('-s', '--species-list', dest='species_list')
    op.add_argument('-i', '--ids-list', dest='ids_list')

    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)

    op.usage = '''Extends an orthogroup database and orthogroups files.
    First argument is a fath to existed Scenario 1 output.

    usage: %s <directory> [--assemblies dir]
                        [--proteomes dir] [--annotations dir] [--ids-list file]
                        [--species-list file]
                        [-t num] [--start-from step]

    First argument <directory> is a fath to existed Scenario 1 output.

    Optional arguments:
    -a --assemblies:     Directory with assemblies in fasta format.

    -g --gbs:            Directory with gb files.

    -p --proteomes:      Directory with fasta (or faa) protein files,
                         named by their reference ids (i.e. NC_005816.1.fasta).
                         Can contain annotations from Prodigal.

    -s --species-list:   File with a list of organism names as in Genbank.

    -i --ids-list:       File with reference ids (will be fetched from Genbank).

    --prot-id-field:     When specifying proteomes, use this fasta id field number
                         to retrieve protein ids (default if 1, like
                         >NC_005816.1|NP_995567.1 ...).
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

    if p.assemblies:
        check_dir(p.species_list)
        p.assemblies = abspath(p.assemblies)

    if p.proteomes:
        check_dir(p.proteomes)
        p.proteomes = abspath(p.proteomes)

    if not isdir(p.directory):
        arg_parse_error('Directory %s does not exist.' % p.directory)

    return p


assemblies_dir = 'assemblies'


def step_finding_genes(assembly_names):
    def run():
        if assembly_names:
            return cmdline('prodigal',
                parameters=[
                    '-i', filtered_assemblies,
                    '-a', predicted_proteomes]),
                   #'-o', join(steps.intermediate_dir, assembly_name + '.gff')]
        else:
            log.info('Skipping.')
            return 0

    filtered_assemblies = \
        [join(steps.intermediate_dir, asm_name + '.fna')
         for asm_name in assembly_names]

    predicted_proteomes = \
        [join(steps.proteomes_dir, asm_name + '.fasta')
         for asm_name in assembly_names]

    return Step(
        'Finding genes',
         run=run,
         req_files=[steps.proteomes_dir, filtered_assemblies],
         prod_files=[predicted_proteomes])


def step_blast_singletones():
    def blast_singletones(singletones_file, new_proteomes_dir):
        log.error('Not implemented: in progress.')
        return 1
        with open(singletones_file) as f:
            for line in f:
                strain, prot_id, gene_id, locus, product, description = line.split()
                # TODO: Blast!
        return 1

    return Step(
       'Blasting singletones',
        run=lambda: blast_singletones(steps.assembly_singletones, new_proteomes_dir),
        req_files=[steps.assembly_singletones])


new_proteomes_dir = join('new_proteomes')


def step_prepare_input(p):
    if p.assemblies:
        assemblies = [
            join(p.assemblies, f)
            for f in listdir(p.assemblies)
            if f and f[0] != '.']
        assembly_names = [
            splitext(basename(asm))[0]
            for asm in assemblies]
        filtered_assemblies = [
            join(assemblies_dir, asm_name + '.fna')
            for asm_name in assembly_names]
        new_proteomes = [
            join(steps.proteomes_dir, asm_name + '.fasta')
            for asm_name in assembly_names]

        def run():
            if not isdir(assemblies_dir): mkdir(assemblies_dir)

            total_successful_filters = 0
            for assembly, filtered_asm in zip(assemblies, filtered_assemblies):
                if filter_assembly(assembly,
                                   filtered_asm,
                                   skip=(4, 7, 10, 23, 32, 38),
                                   skip_after=51) == 0:
                    total_successful_filters += 1
            if total_successful_filters == 0:
                log.error('No correct assemblies.')
                return 1

            for asm, prot, asm_name in zip(
                    filtered_assemblies, new_proteomes, assembly_names):
                res = cmdline('prodigal',
                    parameters=[
                       '-i', asm,
                       '-o', join(steps.intermediate_dir, asm_name),
                       '-a', prot])()
                if res != 0:
                    return res
                log.info('')

            res = adjust_proteomes(new_proteomes, steps.proteomes_dir,
                                   prot_id_field=0)
            if res != 0:
                return res

            if exists(new_proteomes_dir):
                rmtree(new_proteomes_dir)
            if not isdir(new_proteomes_dir):
                mkdir(new_proteomes_dir)
            for prot in new_proteomes:
                copy(prot, join(new_proteomes_dir, basename(prot)))

            return 0

        return Step(
           'Preparing input',
            req_files=assemblies,
            prod_files=new_proteomes,
            run=run)

    elif p.proteomes:
        proteomes = [
            join(p.proteomes, prot)
            for prot in listdir(p.proteomes)
            if prot and prot[0] != '.']

        new_proteomes = [
            join(new_proteomes_dir, basename(prot) + '.fasta')
            for prot in proteomes]

        def run():
            if not isdir(new_proteomes_dir):
                mkdir(new_proteomes_dir)

            for prot in proteomes:
                copy(prot, join(steps.proteomes_dir, basename(prot)))

            for prot in proteomes:
                copy(prot, join(new_proteomes_dir, basename(prot)))

            return 0

        return Step(
           'Preparing input',
            req_files=proteomes,
            prod_files=new_proteomes,
            run=run)


new_good_proteomes = join(steps.intermediate_dir, 'new_good_proteins.fasta')
new_bad_proteomes = join(steps.intermediate_dir, 'new_bad_proteins.fasta')


def filter_new_proteomes(new_proteomes_dir, min_length=10, max_percent_stop=20):
    return Step(
        'Filtering new proteomes',
         run=cmdline(join(steps.orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
             parameters=[
                 new_proteomes_dir,
                 min_length, max_percent_stop,
                 new_good_proteomes,
                 new_bad_proteomes]),
         req_files=[new_proteomes_dir],
         prod_files=[
             new_good_proteomes,
             new_bad_proteomes])


def main(args):
    register_ctrl_c()

    p = parse_args(args)
    log_fpath = set_up_logging(p.debug, p.directory)
    log.info('python ' + basename(__file__) + ' ' + ' '.join(args))
    log.info('')

    check_and_install_tools(p.debug, log_fpath)
    set_up_config()

    start_from, start_after = get_starting_step(p.start_from, join(p.directory, log_fname))

    working_dir = p.directory
    log.info('Changing to %s' % working_dir)
    chdir(working_dir)

    if not exists('intermediate'):
        arg_parse_error('You need to run Scenario 1 on this directory first.')

    workflow = Workflow(working_dir, id=make_workflow_id(working_dir))
    log.info('Workflow id is "' + workflow.id + '"')
    log.info('')
    suffix = '_' + workflow.id

    workflow.extend([
        step_prepare_input(p),
        steps.filter_proteomes(
            min_length=int(p.min_length),
            max_percent_stop=int(p.max_percent_stop)),
        filter_new_proteomes(
            new_proteomes_dir,
            min_length=int(p.min_length),
            max_percent_stop=int(p.max_percent_stop)),
        steps.make_blast_db(),
        steps.blast(p.threads, new_good_proteomes, evalue=float(p.evalue)),
        steps.parse_blast_results(),
        steps.clean_database(suffix),
        steps.install_schema(suffix),
        steps.load_blast_results(suffix),
        steps.find_pairs(suffix),
        steps.dump_pairs_to_files(suffix),
        steps.mcl(),
        steps.step_save_orthogroups(new_proteomes_dir),
        step_blast_singletones(),
    ])

    result = workflow.run(start_after, start_from,
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
    main(sys.argv[1:])