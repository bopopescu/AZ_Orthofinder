#!/usr/bin/env python

from random import randint
from shutil import rmtree, copyfile, copy
from subprocess import call
from os import chdir, getcwd, listdir, mkdir, remove
from os.path import join, realpath, isfile, exists, dirname, basename, normpath, isdir
from itertools import ifilterfalse, ifilter, izip, count
from mysql.connector import errorcode
import mysql.connector
import sys
import logging
from az_orthofinder.src.Workflow import Step, Workflow

import utils
from src.fetch_annotations import fetch_annotations_ids, fetch_annotations_species_name_ftp
from src.make_proteomes import make_proteomes
from src.db_connection import DbCursor

log = utils.log
orthomcl_config = utils.orthomcl_config
orthomcl_bin_dir = utils.orthomcl_bin_dir


def parse_args(args):
    import argparse
    op = argparse.ArgumentParser(description='Find ortholougs for an assembly '
                                             'based on existing database.')

    op.add_argument('-a', '--assembly', dest='assembly', help='Fasta file with contigs.')

    op.add_argument('-b', '--blast-results', dest='blast_results',
                    help='Tab-separated blast file.')

    op.add_argument('-w', '--workflow-id', dest='workflow_id', default=None)

    op.add_argument('-o', dest='out_dir',
                    help='The directory that will contain the resulting group.txt file, '
                         'as well as intermediate results.')

    op.add_argument('-w', '--overwrite', dest='overwrite', action='store_true', default=False,
                    help='By default, the tool reuses existing intermediate results.'
                         'This option makes the tool preliminarily remove existing data')

    op.add_argument('--start-from', dest='start_from', default=0,
                    help='Start from this step. Either name (see log.txt) or "uselog".'
                         'If "uselog", find the last "Done" record in log.txt file.')

    op.add_argument('-t', '--threads', dest='threads', default=1,
                    help='Number of threads to run Blast.')

    op.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)

    op.add_argument('--ask-each-step',
                    dest='ask_each_step', action='store_true', default=False,
                    help='Ask user every time before proceed to next step.')

    op.add_argument('--proxy', dest='proxy', default=None, help='Set up proxy for FTP.')

    params = op.parse_args(args)
    if params.start_from == 'uselog':
        if not isfile(join(params.out_dir, utils.log_fname)):
            print >> sys.stderr, 'No %s in %s. Either check your path, or ' \
                                 'change the --start-from=uselog option' %\
                                 (utils.log_fname, params.out_dir)
            exit(1)

    if params.workflow_id is None:
        params.workflow_id = basename(params.out_dir)[:4]

    if params.species_file is None and params.ref_ids is None:
        print >> sys.stderr, 'Either --ref-ids or --species-file should be specified.'
        exit(1)

    return params



def run_workflow(working_dir,
                 workflow_id,
                 old_blast_out,
                 overwrite,
                 ask_before=False,
                 start_after=None,
                 start_from=None,
                 threads=1,
                 proxy=None,
                 inflation=1.5,
                 first_id=1000):

    if not exists('intermediate'):
        mkdir('intermediate')

    proteomes_dir       = 'proteomes'
    annotations_dir     = 'annotations'
    sql_log             = 'intermediate/log.sql'
    good_proteins       = 'intermediate/good_proteins.fasta'
    poor_proteins       = 'intermediate/poor_proteins.fasta'
    new_blast_out       = 'intermediate/blasted.tsv'
    similar_sequences   = 'intermediate/similar_sequences.txt'
    pairs_log           = 'intermediate/orthomclpairs.log'
    mcl_input           = 'mclInput'
    mcl_output          = 'mclOutput'
    groups_file         = 'groups.txt'
    singletons_file     = 'singletons.txt'

    with open(orthomcl_config) as f:
        conf = dict(l.split('=', 1) for l in f.readlines() if l[0] != '#')
        ortholog_table = conf['orthologTable'].strip() + '_' + workflow_id
        in_paralog_table = conf['inParalogTable'].strip() + '_' + workflow_id
        coortholog_table = conf['coOrthologTable'].strip() + '_' + workflow_id
        similar_sequeces_table = conf['similarSequencesTable'].strip() + '_' + workflow_id
        inter_taxon_match_view = conf['interTaxonMatchView'].strip() + '_' + workflow_id
        best_hit_table = 'BestHit' + '_' + workflow_id
        best_hit_taxon_score_table = 'BestQueryTaxonScore' + '_' + workflow_id

    steps = [
        Step('Making proteins',
             cmd=make_proteomes,
             req_files=[annotations_dir],
             prod_files=[proteomes_dir],
             parameters=[
                annotations_dir,
                workflow_id,
                proteomes_dir,
                ]),

        Step('Filtering fasta',
             cmd=join(orthomcl_bin_dir, 'orthomclFilterFasta.pl'),
             req_files=[proteomes_dir],
             prod_files=[good_proteins, poor_proteins],
             parameters=[
                proteomes_dir,
                10,
                20,
                good_proteins,
                poor_proteins,
                ]),

        Step('Blasting all vs. all',
             cmd='blastp',
             req_files=[proteomes_dir, good_proteins],
             prod_files=[blast_out],
             parameters=[
                '-query', good_proteins,
                '-subject', good_proteins,
                '-out', blast_out,
                '-outfmt', 6,  # tabular
                '-evalue', 1e-5,
                '-num_descriptions', 10000,  # don't care value
                '-num_alignments', 10000,  # don't care value,
                '-num_threads', threads,
                '-dbsize', blast_dbsize,
                ]),

        Step('Parsing blast results',
             cmd=join(orthomcl_bin_dir, 'orthomclBlastParser.pl'),
             req_files=[proteomes_dir, blast_out],
             prod_files=[similar_sequences],
             parameters=[
                blast_out,
                proteomes_dir],
             stdout=similar_sequences),

        Step('Installing schema',
             cmd=join(orthomcl_bin_dir, 'orthomclInstallSchema.pl'),
             req_files=[orthomcl_config],
             prod_files=[],  # drops and creates SimilarSequences, InParalog,
                             # Ortholog, CoOrtholog, InterTaxonMatch view
             prod_tables=[
                ortholog_table,
                in_paralog_table,
                coortholog_table,
                similar_sequeces_table,
                inter_taxon_match_view],
             parameters=[
                orthomcl_config,
                sql_log,
                ('_' + workflow_id if workflow_id else ''),
                ]),

        Step('Loading blast results into the database',
             cmd=join(orthomcl_bin_dir, 'orthomclLoadBlast.pl'),
             req_files=[orthomcl_config, similar_sequences],  # and initialized database
             prod_files=[],  # loads blast results into the db
             parameters=[
                orthomcl_config,
                similar_sequences,
                ('_' + workflow_id) if workflow_id else ''
                ]),

        Step('Finding pairs',
             cmd=join(orthomcl_bin_dir, 'orthomclPairs.pl'),
             req_files=[orthomcl_config],  # and initialized database
             prod_files=[],  # populates InParalog, Ortholog, CoOrtholog
             parameters=[
                orthomcl_config,
                pairs_log,
                'cleanup=no',
                #'startAfter=useLog' if not overwrite else '',
                ('suffix=_' + workflow_id) if workflow_id else '',
                ]),

        Step('Dump pairs files',
             cmd=join(orthomcl_bin_dir, 'orthomclDumpPairsFiles.pl'),
             req_files=[orthomcl_config],  # and populated InParalog, Ortholog, CoOrtholog tables
             prod_files=[
                mcl_input,
                'pairs',
                'pairs/potentialOrthologs.txt',
                'pairs/potentialInparalogs.txt',
                'pairs/potentialCoorthologs.txt'],
             parameters=[
                orthomcl_config,
                working_dir,  # current dir. no effect, just to use next positional argument
                ('_' + workflow_id) if workflow_id else '']),

        Step('MCL',
             cmd='mcl',
             req_files=[mcl_input],
             prod_files=[mcl_output],
             parameters=[
                mcl_input,
                '--abc',
                '-I', str(inflation),
                '-o', mcl_output,
                ]),

        Step('MCL groups to files',
             cmd=join(orthomcl_bin_dir, 'orthomclMclToGroups.pl'),
             req_files=[mcl_output],
             prod_files=[groups_file],
             parameters=[
                workflow_id + '_',
                str(first_id)],
             stdin='mclOutput',
             stdout=groups_file),

        Step('MCL singletones to files',
             cmd=join(orthomcl_bin_dir, 'orthomclSingletons.pl'),
             req_files=[good_proteins, groups_file],
             prod_files=[singletons_file],
             parameters=[
                good_proteins,
                groups_file],
             stdout=singletons_file),
    ]

    workflow = Workflow(steps)
    return workflow.run(start_after, start_from, overwrite, ask_before)


def main(args):
    params = parse_args(args)

    if not isdir(params.out_dir):
        mkdir(params.out_dir)

    copy(params.blast_results, params.out_dir)

    start_after = None
    start_from = params.start_from
    if params.start_from == 'uselog':
        with open(join(params.out_dir, utils.log_fname)) as f:
            for l in f.readlines():
                if len(l) > 41 and l[36:40] == 'Done':
                    start_after = l[41:].strip()
                    start_from = None

    utils.set_up_logging(params.debug, params.out_dir)

    log.info('Changing to %s' % params.out_dir)
    chdir(params.out_dir)

    run_workflow(working_dir=params.out_dir,
                 workflow_id=params.workflow_id,
                 old_blast_out=params.blast_results,
                 overwrite=params.overwrite or start_after != 0 or start_from != 0,
                 ask_before=params.ask_each_step,
                 start_after=start_after,
                 start_from=start_from,
                 threads=params.threads,
                 proxy=params.proxy)


if __name__ == '__main__':
    main(sys.argv[1:])