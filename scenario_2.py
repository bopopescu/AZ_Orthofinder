#!/usr/bin/env python
from collections import namedtuple
from genericpath import isfile
from shutil import copyfile, rmtree, copy, copytree

import sys
import logging
from os import chdir, mkdir, getcwd, listdir, symlink, makedirs, rmdir, remove, environ
from os.path import join, exists, isdir, dirname, realpath,\
    basename, splitext, abspath, expanduser
import urllib2
from Bio import SeqIO
from src.fetch_annotations import fetch_annotations_for_ids
from src.process_assembly import filter_assembly
from src.make_proteomes import adjust_proteomes, make_proteomes
from src import steps
import time

from src.utils import which, make_workflow_id, read_list, \
    set_up_config, get_starting_step, register_ctrl_c, \
    check_installed_tools, test_entrez_conn, \
    check_install_mcl, check_perl_modules, check_and_install_tools, test_blast_conn
from src.parse_args import arg_parse_error, check_file,\
    check_dir, add_common_arguments, check_common_args
from src.logger import set_up_logging, add_file_handler
from src.Workflow import Workflow, Step, cmdline

from src import config

from src.config import log_fname
log = logging.getLogger(log_fname)

script_path = dirname(realpath(__file__))

max_attempts = 1


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find groups of orthologous genes.')

    #op.add_argument(dest='directory', required=False)
    op.add_argument('-s1o', dest='directory', required=True)
    op.add_argument('-s2o', '-o', dest='out_dir', required=False)

    op.add_argument('-a', '--assemblies', dest='assemblies')
    op.add_argument('-g', '--annotations', '--gbs', dest='annotations')
    op.add_argument('-p', '--proteomes', '--proteins', dest='proteomes')
    op.add_argument('-i', '--ids', '--ids-list', dest='ids_list')

    op.add_argument('--prot-id-field', dest='prot_id_field', default=1)
    op.add_argument('--skip-blast-singletones', dest='blast_singletones',
                    action='store_false', default=True)
    op.add_argument('--blastdb', '--blast-db', dest='blastdb')

    #-o:                  Output directory (if not specified, the input directory will be used).

    op.usage = '''Extends an orthogroup database and orthogroups files.
First argument is a fath to existed Scenario 1 output.

Test runs:
    python scenario_2.py -s1o test_ids -s2o test_ids_new_ids --ids test_input/new_ids.txt

    python scenario_2.py -s1o test_proteomes -s2o test_prots_new_prots --proteomes test_input/new_proteins

Usage: %s -s1o <scenario_1 directory> -s2o <out_dir> [--assemblies dir] [--proteomes dir]
                                 [--gbs dir] [--ids-list file] [--species-list file]
                                 [-t num] [--start-from step] [--blast-singletones] [--blast-db]

    -s1o                 Path to existed Scenario 1 output.

    -s2o                 Output directory (optional, if ommited, the input directory will be used).

Optional arguments:

    -a --assemblies:     Directory with assemblies in fasta format.

    -g --gbs:            Directory with gb files.

    -p --proteomes:      Directory with fasta (or faa) protein files,
                         named by their reference ids (i.e. NC_005816.1.fasta).
                         Can contain annotations from Prodigal.

    -i --ids-list:       File with reference ids (will be fetched from Genbank).

    --prot-id-field:     When specifying proteomes, use this fasta id field number
                         to retrieve protein ids (default if 1, like
                         >NC_005816.1|NP_995567.1 ...).

    --blast-singletones  Search newly added proteins agains NCBI database, if they did
                         not fit any group with known proteins.

    --blast-db           Local Blast database path. If not set, remote NCBI will be used.
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
        check_dir(expanduser(p.assemblies))
        p.assemblies = abspath(expanduser(p.assemblies))

    if p.proteomes:
        check_dir(expanduser(p.proteomes))
        p.proteomes = abspath(expanduser(p.proteomes))

    if p.ids_list:
        check_file(expanduser(p.ids_list))
        p.ids_list = abspath(expanduser(p.ids_list))

    if p.annotations:
        check_dir(expanduser(p.annotations))
        p.annotations = abspath(expanduser(p.annotations))

    if p.blastdb:
        p.blastdb = abspath(expanduser(p.blastdb))

    if not isdir(expanduser(p.directory)):
        arg_parse_error('Directory %s does not exist.' % p.directory)
    p.directory = abspath(expanduser(p.directory))

    if not p.out_dir:
        p.out_dir = p.directory
        #arg_parse_error('Specify output directory with -o.')
    if isfile(expanduser(p.out_dir)):
        arg_parse_error('%s is a file' % p.out_dir)
    p.out_dir = abspath(expanduser(p.out_dir))

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
        [join(config.intermediate_dir, asm_name + '.fna')
         for asm_name in assembly_names]

    predicted_proteomes = \
        [join(config.proteomes_dir, asm_name + '.fasta')
         for asm_name in assembly_names]

    return Step(
        'Finding genes',
         run=run,
         req_files=[config.proteomes_dir, filtered_assemblies],
         prod_files=[predicted_proteomes])


blasted_singletones_dir = 'blasted_singletones'


def blast_online(rec, result_xml_fpath):
    from Bio.Blast import NCBIWWW

    retrying = False
    to_the_next = False
    attempt_number = 1
    while True:
        try:
            print rec.format('fasta')
            result_xml_f = NCBIWWW.qblast('blastp', 'refseq_protein', rec.format('fasta'),
                                        hitlist_size=10)
            with open(result_xml_fpath, 'w') as save_f:
                save_f.write(result_xml_f.read())

        except urllib2.HTTPError as e:
            log.warn('     Warning: could not blast through web. HTTPError: %s. Code %s. '
                     '(You can press Ctrl+C to interrupt and continue later).'
                     % (e.msg, str(e.code)))
            retrying = True

        except urllib2.URLError, e:
            log.warn('     Warning: could not blast through web. URLError: %s. '
                     '(You can press Ctrl+C to interrupt and continue later).'
                     % (e.args))
            retrying = True

        except (KeyboardInterrupt, SystemExit, GeneratorExit):
            if retrying:
                log.info('     If you restart from this step and do not remove the "%s" directory, '
                         'the process will continue from here.' % blasted_singletones_dir)
            return 1
        time.sleep(2)

        if retrying:
            if attempt_number >= max_attempts:
                to_the_next = True
                break
            else:
                attempt_number += 1
                log.info('     Attempt %d/3' % attempt_number)

    if not to_the_next:
        time.sleep(0.5)
    return 0


def search_for_best_hit(rec, res_xml_fpath, short_fpath):
    from Bio.Blast import NCBIXML

    LEN_FACTOR = 0.05
    BIG_EVALUE = 2
    BestHits = namedtuple('BestHits', 'score, evalue, alignments, hits')
    best_hits = BestHits(0, BIG_EVALUE, set(), set())

    with open(res_xml_fpath) as full_xml_f, \
         open(short_fpath, 'w') as short_f:
        short_f.write(rec.description + '\n')
        short_f.write(str(rec.seq) + '\n\n')

        blast_record = next(NCBIXML.parse(full_xml_f))

        for i, alignment in enumerate(blast_record.alignments):
            short_f.write(str(i + 1) + '. Alignment\n'
                          '   Title: ' + alignment.title + '\n'
                          '   Length: ' + str(alignment.length) + '\n'
                          '   Accession: ' + alignment.hit_id + '\n')

            #if 'protein' in alignment.title:
            for hsp in alignment.hsps:
                short_f.write(
                    '     Hit score: ' + str(hsp.score) + '\n'
                    '     Hit expect value: ' + str(hsp.expect) + '\n'
                    '     Hit query (starts at ' + str(hsp.query_start) + '):\n     ' + hsp.query + '\n'
                    '     Hit match:\n     ' + hsp.match + '\n'
                    '     Hit subject (starts at ' + str(hsp.sbjct_start) + ':\n     ' + hsp.sbjct + '\n\n')

                if hsp.expect != 0:
                    if hsp.expect == best_hits.evalue:
                        best_hits.hits += hsp
                        best_hits.alignments += alignment
                    if hsp.expect < best_hits.evalue:
                        best_hits = BestHits(hsp.score, hsp.expect, {alignment}, {hsp})
                else:
                    if (len(rec.seq) / len(hsp.match)) - 1 > LEN_FACTOR:
                        log.debug('     Evaue is 0 and lengths do not match: '
                                  'len(rec.seq)/len(match) - 1 = %f, witch is greater than '
                                  'the threshold of %f; Title = %s' %
                                  ((len(rec.seq) / len(hsp.match)) - 1 , LEN_FACTOR, alignment.title))
                    else:
                        log.debug('     len(rec.seq)/len(match) - 1 = ' +
                                  str((len(rec.seq) / len(hsp.match)) - 1))
                        if hsp.score == best_hits.score:
                            best_hits.hits.add(hsp)
                            best_hits.alignments.add(alignment)
                        if hsp.score > best_hits.score:
                            best_hits = BestHits(hsp.score, hsp.expect, {alignment}, {hsp})

        if best_hits.hits:
            log.info('     e-value: ' + str(best_hits.evalue))
            log.info('     score:   ' + str(best_hits.score))
            for hit, alignment in zip(best_hits.hits, best_hits.alignments):
                log.info('     hits:')
                log.info('       title:     ' + alignment.title)
                log.info('       accession: ' + alignment.hit_id)
                log.info('       length:    ' + str(alignment.length))
                log.info('       ' + hit.query[:75] + '...')
                log.info('       ' + hit.match[:75] + '...')
                log.info('       ' + hit.sbjct[:75] + '...')
        else:
            log.warning('     No hits for ' + rec.id)
    log.info('     Saved to ' + short_fpath)
    log.info('')


def process_record(rec, group_i, blastdb, threads):
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastpCommandline

    log.info('     Reading sequence ' + rec.id)

    # Blasting against NCBI
    res_xml_fpath = join(
        blasted_singletones_dir,
        'group_' + str(group_i + 1) + '_refseq_blasted_' + rec.id.replace('|', '__') + '.xml')
    short_fpath = join(
        blasted_singletones_dir,
        'group_' + str(group_i + 1) + '_refseq_blasted_' + rec.id.replace('|', '__') + '_summary.txt')

    do_blast = True
    if isfile(res_xml_fpath):
        with open(res_xml_fpath) as full_xml_f:
            try:
                next(NCBIXML.parse(full_xml_f))
            except:
                log.info(
                    '     File %s is empty (or could no be read). '
                    'The sequences is going to be blasted again. the file.' % abspath(res_xml_fpath))
                remove(res_xml_fpath)
            else:
                log.info('     Reading results from ' + abspath(res_xml_fpath))
                do_blast = False

    if do_blast:
        log.info('     Blasting against the refseq_proteins database.')
        log.info('     Writing result to ' + abspath(res_xml_fpath))

        if blastdb:
            blast_cmdline = NcbiblastpCommandline(
                db=blastdb,
                outfmt=5,
                out=res_xml_fpath,
                num_threads=threads)
            out, err = blast_cmdline(stdin=str(rec.seq))
        else:
            res = blast_online(rec, res_xml_fpath)
            if res == 1:
                return 1

    search_for_best_hit(rec, res_xml_fpath, short_fpath)


def step_blast_singletones(threads, blast_singletones=True, blastdb=None, debug=False, rewrite=False):
    def run(singletones_file, new_proteomes_dir):
        #environ["http_proxy"] = "http://192.168.0.2:3128"

        #if not blast_singletones:
        #    raw_input('Blast added proteomes? ' +
        #              (('Local database ' + blastdb + ' will be used')
        #               if blastdb
        #               else 'Remote database will be used.')

        if isdir(config.singletone_dir) and \
            [join(config.singletone_dir, fname)
             for fname in listdir(config.singletone_dir)
             if fname and fname[0] != '.']:
            pass
        else:
            log.error('   No singletones in additional genomes, skipping the step.')
            log.debug('   ' + realpath(config.singletone_dir))
            return 0

        if blastdb:
            log.info('   Using local NCBI database: ' + blastdb)
        else:
            #if test_blast_conn():
            if not blast_singletones:
                try:
                    raw_input('Blast added proteomes? '
                              'Remote database will be used. '
                              'Press any key to overwrite and continue, '
                              'or Ctrl+C to interrupt.\n> ')
                except (EOFError, KeyboardInterrupt, SystemExit, GeneratorExit):
                    exit(1)
            log.info('   Using remote NCBI database.')
            #else:
            #    log.error('   No Blast database and no Internet connection '
            #              'to use the remote NCBI database. Please, provide '
            #              'a path to blast database (with the --blast-db option) '
            #              'or verify your Internet connection.')
            #    return 1

        if rewrite and exists(blasted_singletones_dir):
            rmtree(blasted_singletones_dir)
        if not isdir(blasted_singletones_dir):
            mkdir(blasted_singletones_dir)

        for group_i, group_singletones_file in enumerate(
                (join(config.singletone_dir, fname)
                 for fname in listdir(config.singletone_dir)
                 if fname and fname[0] != '.')):
            log.debug('   Group ' + str(group_i + 1) + '. ' + group_singletones_file)

            for rec in SeqIO.parse(group_singletones_file, 'fasta'):
                res = process_record(rec, group_i, blastdb, threads)
                if res == 1:
                    return 1

        return 0

    return Step(
        'Blasting singletones',
        run=lambda: run(config.assembly_singletones_file, new_proteomes_dir))


new_proteomes_dir = 'new_proteomes'
new_annotations_dir = 'new_annotations'


def filter_dublicated_proteomes(prot_dir, new_files):
    prot_names = \
        [splitext(prot)[0] for prot in listdir(prot_dir) if prot and prot[0] != '.']
    new_prot_names = \
        [splitext(basename(prot))[0] for prot in new_files if prot and prot[0] != '.']

    filtered_new_prots = []
    for new_prot_name, new_prot_path in zip(new_prot_names, new_files):
        if new_prot_name in prot_names:
            log.warn('   Proteome "%s" is already considered (check "%s").' \
                      % (new_prot_name, realpath(prot_dir)))
        else:
            filtered_new_prots.append(new_prot_path)

    return filtered_new_prots


all_considered_warning = '   Notice: all proteomes are already considered in this directory.\n' + \
                         'If you are sure the input is different from the proteomes in the %s directory, ' + \
                         'you will need to rename the input files.\nNote that proteomes are identified ' + \
                         'by their filenames, so it is not desired to just removed those files and restart.'

def step_prepare_input(p):
    if p.assemblies:
        def run():
            assemblies = [
                join(p.assemblies, f)
                for f in listdir(p.assemblies)
                if f and f[0] != '.']

            if isdir(config.proteomes_dir):
                assemblies = filter_dublicated_proteomes(config.proteomes_dir, assemblies)
                if assemblies == []:
                    log.warn(all_considered_warning % config.proteomes_dir)
                    exit(1)

            assembly_names = [
                splitext(basename(asm))[0]
                for asm in assemblies]
            filtered_assemblies = [
                join(assemblies_dir, asm_name + '.fna')
                for asm_name in assembly_names]
            new_proteomes = [
                join(config.proteomes_dir, asm_name + '.fasta')
                for asm_name in assembly_names]

            if not isdir(assemblies_dir): mkdir(assemblies_dir)
            log.debug('   Created assemblies_dir ' + assemblies_dir)

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
                       '-o', join(config.intermediate_dir, asm_name),
                       '-a', prot])()
                if res != 0:
                    return res
                log.info('')

            res = adjust_proteomes(new_proteomes, config.proteomes_dir,
                                   prot_id_field=0)
            if res != 0:
                return res

            # Recreate new_proteomes_directory
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

            if isdir(config.proteomes_dir):
                input_proteomes = filter_dublicated_proteomes(config.proteomes_dir, input_proteomes)
                if input_proteomes == []:
                    log.warn(all_considered_warning % config.proteomes_dir)
                    exit(1)

            new_prot_names = [splitext(basename(prot))[0] for prot in input_proteomes]

            new_proteomes_dir = 'new_proteomes'
            new_annotations_dir = 'new_annotations'
            new_proteomes = [
                join(new_proteomes_dir, prot_name + '.fasta')
                for prot_name in new_prot_names]

            new_proteomes_dir = join(getcwd(), new_proteomes_dir)
            if not isdir(new_proteomes_dir):
                mkdir(new_proteomes_dir)

            for prot_from, prot_to in zip(input_proteomes, new_proteomes):
                copy(prot_from, prot_to)
                copy(prot_from, join(config.proteomes_dir, basename(prot_to)))

            return 0

        return Step(
           'Preparing input',
            run=run)

    elif p.ids_list:
        def run():
            if not test_entrez_conn():
                log.error('No internet connection: cannot fetch annotations.')
                return 4

            log.debug('   Using ref ids: ' + str(p.ids_list))
            ref_ids = read_list(p.ids_list)
            res = fetch_annotations_for_ids(new_annotations_dir, ref_ids)
            if res != 0: return res
            res = make_proteomes(new_annotations_dir, new_proteomes_dir)
            if res != 0: return res
            for fname in listdir(new_annotations_dir):
                if fname[0] != '.':
                    copy(join(new_annotations_dir, fname), config.annotations_dir)
            for fname in listdir(new_proteomes_dir):
                if fname[0] != '.':
                    copy(join(new_proteomes_dir, fname), config.proteomes_dir)

            return 0

        return Step(
           'Preparing input',
            run=run)


new_good_proteomes = join(config.intermediate_dir, 'new_good_proteins.fasta')
new_bad_proteomes = join(config.intermediate_dir, 'new_bad_proteins.fasta')


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

    try:
        if not exists(join(p.directory, 'intermediate')):
            arg_parse_error('You need to run Scenario 1 on this directory first.')

        if not p.out_dir:
            p.out_dir = p.directory

        working_dir = p.out_dir

        if p.out_dir != p.directory:
            log_text = ''

            if isdir(p.out_dir):
                if not p.overwrite:
                    files = [f for f in listdir(p.out_dir) if f and f[0] != '.']
                    #log.debug(files)
                    if files:
                        print('The output directory exists. Do you want to overwrite it? ' +
                              '(You can run with the --overwrite option to avoid this warning.)')
                        try:
                            raw_input('Press any key to overwrite and continue, or Ctrl+C to interrupt.\n> ')
                        except (EOFError, KeyboardInterrupt, SystemExit, GeneratorExit):
                            return 1
                if exists(join(p.out_dir, log_fname)):
                    with open(join(p.out_dir, log_fname)) as log_f:
                        log_text = log_f.read()
                rmtree(p.out_dir)

            makedirs(p.out_dir)
            rmdir(p.out_dir)
            copytree(p.directory, p.out_dir)
            if isfile(join(p.out_dir, log_fname)):
                remove(join(p.out_dir, log_fname))
            chdir(p.out_dir)
            if log_text:
                with open(join(p.out_dir, log_fname), 'w') as log_f:
                    log_f.write(log_text)

        log_fpath = set_up_logging(p.debug, p.out_dir, 'a')
        log.info('python ' + basename(__file__) + ' ' + ' '.join(args))
        log.info('')
        check_and_install_tools(p.debug, log_fpath)
        set_up_config(working_dir)

        start_from, start_after = get_starting_step(p.start_from, join(working_dir, log_fname))

        log.info('Changing to %s' % working_dir)
        if not isdir(working_dir):
            makedirs(working_dir)
        chdir(working_dir)

        workflow = Workflow(working_dir, id=make_workflow_id(working_dir),
                            cmdline_args=['python', basename(__file__)] + args)
        log.info('Workflow id is "' + workflow.id + '"')
        log.info('')

        with open(config.config_file) as f:
            conf = dict(l.strip().lower().split('=', 1)
                        for l in f.readlines() if l.strip() and l.strip()[0] != '#')
            if conf['db_vendor'] == 'sqlite':
                suffix = ''
            else:
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
            steps.step_save_orthogroups(new_proteomes_dir if not p.ids_list and p.blast_singletones else None)
        ])
        if not p.ids_list:
            workflow.extend([step_blast_singletones(p.threads, p.blast_singletones, p.blastdb, p.debug, p.overwrite)])

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
                log.info('Groups in short format are in ' + join(working_dir, config.short_orthogroups_file))

        return result

    except (KeyboardInterrupt, SystemExit, GeneratorExit):
        return 1

    except Exception as e:
        log.error('')
        log.exception('Unexpected error!')
        return 2


if __name__ == '__main__':
    main(sys.argv[1:])