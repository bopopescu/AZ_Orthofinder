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
    op.add_argument('-p', '--proteomes', '--proteins', dest='proteomes')
    op.add_argument('-s', '--species-list', dest='species_list')
    op.add_argument('-i', '--ids-list', dest='ids_list')

    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)
    op.add_argument('--blastdb', dest='blastdb')

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

    --blast-db           Local Blast database path. If not set, NCBI will be used.
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


blasted_singletones_dir = 'blasted_singletones'


def step_blast_singletones(blastdb=None, debug=False):
    def run(singletones_file, new_proteomes_dir):
        from Bio.Blast import NCBIWWW, NCBIXML
        from Bio import SeqIO

        if not blastdb:
            if test_internet_conn():
                log.info('   Using remote NCBI database.')
            else:
                log.error('   No Blast database and no internet connection'
                          ' to use the remote NCBI database.')
                return 1

            #if exists(blasted_singletones_dir):
            #    rmtree(blasted_singletones_dir)
            #if not isdir(blasted_singletones_dir):
            #    mkdir(blasted_singletones_dir)

            for group_singletones_file in (join(steps.singletone_dir, fname)
                                           for fname in listdir(steps.singletone_dir)
                                           if fname and fname[0] != '.'):
                log.debug('   ' + group_singletones_file)
                rec = next(SeqIO.parse(group_singletones_file, 'fasta'))
                log.debug('     ' + rec.id)

                # Blasting against NCBI
                save_fpath = join(blasted_singletones_dir, 'refseq_blasted_' + rec.id + '.xml')
                if debug and isfile(save_fpath):
                    pass
                else:
                    result_handle = NCBIWWW.qblast('blastp', 'refseq_protein', rec.format('fasta'))
                    with open(save_fpath, 'w') as save_f:
                        save_f.write(result_handle.read())

                # Searching best hit
                best_alignment, best_hit, best_e = None, None, 2
                with open(save_fpath) as result_handle:
                    blast_record = NCBIXML.read(result_handle)
                    for alignment in blast_record.alignments:
                        if 'protein' in alignment.title:
                            for hsp in alignment.hsps:
                                if hsp.expect < best_e:
                                    best_alignment, best_hit, best_e = alignment, hsp, hsp.expect

                    if best_hit:
                        log.debug('     sequence:' + best_alignment.title)
                        log.debug('     accession: ' + best_alignment.hit_id)
                        log.debug('     length:' + str(best_alignment.length))
                        log.debug('     e value:' + str(best_hit.expect))
                        log.debug('     ' + best_hit.query[:75] + '...')
                        log.debug('     ' + best_hit.match[:75] + '...')
                        log.debug('     ' + best_hit.sbjct[:75] + '...')
                        log.debug('')
                    else:
                        log.warning('     No protein hits for ' + rec.id)

        return 0

    return Step(
       'Blasting singletones',
        run=lambda: run(steps.assembly_singletones_file, new_proteomes_dir),
        req_files=[steps.assembly_singletones_file])


new_proteomes_dir = join('new_proteomes')


def filter_dublicated_proteomes(prot_dir, new_files):
    prot_names = \
        [splitext(prot)[0] for prot in listdir(prot_dir) if prot and prot[0] != '.']
    new_prot_names = \
        [splitext(basename(prot))[0] for prot in new_files if prot and prot[0] != '.']

    filtered_new_prots = []
    for new_prot_name, new_prot_path in zip(new_prot_names, new_files):
        if new_prot_name in prot_names:
            log.warn(('   Proteome "%s" is already considered (check the "' % new_prot_name)
                     + prot_dir + '" directory)')
        else:
            filtered_new_prots.append(new_prot_path)

    return filtered_new_prots


def step_prepare_input(p):
    if p.assemblies:
        def run():
            assemblies = [
                join(p.assemblies, f)
                for f in listdir(p.assemblies)
                if f and f[0] != '.']

            if isdir(steps.proteomes_dir):
                assemblies = filter_dublicated_proteomes(steps.proteomes_dir, assemblies)
                if assemblies == []:
                    log.warn('   Notice: All proteomes are already considered in this directory.'
                             ' If you are sure the input is different from the proteomes in the %s directory, '
                             ' You will need to rename the input files.' % steps.proteomes_dir)
                    exit(1)

            assembly_names = [
                splitext(basename(asm))[0]
                for asm in assemblies]
            filtered_assemblies = [
                join(assemblies_dir, asm_name + '.fna')
                for asm_name in assembly_names]
            new_proteomes = [
                join(steps.proteomes_dir, asm_name + '.fasta')
                for asm_name in assembly_names]

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
            run=run)

    elif p.proteomes:
        def run():
            input_proteomes = [
                join(p.proteomes, prot)
                for prot in listdir(p.proteomes)
                if prot and prot[0] != '.']

            if isdir(steps.proteomes_dir):
                input_proteomes = filter_dublicated_proteomes(steps.proteomes_dir, input_proteomes)
                if input_proteomes == []:
                    log.warn('Notice: All proteomes are already considered in this directory. '
                             'If you are sure the input is different from the proteomes in the %s directory, '
                             'You will need to rename the input files.' % steps.proteomes_dir)
                    exit(1)

            new_prot_names = [splitext(basename(prot))[0] for prot in input_proteomes]

            new_proteomes = [
                join(new_proteomes_dir, prot_name + '.fasta')
                for prot_name in new_prot_names]

            if not isdir(new_proteomes_dir):
                mkdir(new_proteomes_dir)

            for prot_from, prot_to in zip(input_proteomes, new_proteomes):
                copy(prot_from, prot_to)
                copy(prot_from, join(steps.proteomes_dir, basename(prot_to)))

            return 0

        return Step(
           'Preparing input',
            run=run)


new_good_proteomes = join(steps.intermediate_dir, 'new_good_proteins.fasta')
new_bad_proteomes = join(steps.intermediate_dir, 'new_bad_proteins.fasta')


def filter_new_proteomes(new_proteomes_dir, min_length=10, max_percent_stop=20):
    return Step(
       'Filtering new proteomes',
        run=cmdline(
            join(steps.orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
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
        steps.mcl(p.debug),
        steps.step_save_orthogroups(new_proteomes_dir),
        step_blast_singletones(p.blastdb, p.debug),
    ])

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
    main(sys.argv[1:])