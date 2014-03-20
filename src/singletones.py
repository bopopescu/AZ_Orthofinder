from os.path import isdir, join
import urllib2
import logging
from shutil import copyfile, rmtree, copy, copytree
from os import listdir
from os.path import join, exists, isdir, isfile, dirname, realpath,\
    basename, splitext, abspath, expanduser
from os import chdir, mkdir, getcwd, listdir, symlink, makedirs, rmdir, remove
from Bio import SeqIO
import time

from src import config
from src.Workflow import Step, cmdline
from src.config import log_fname, config_file
log = logging.getLogger(log_fname)


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
            if not isfile(blastdb + '.pal'):
                log.error('Incorrect path blast database. Make sure '
                          '%s file exists.' % (blastdb + '.pal'))
                return 1

            os.environ['BLASTDB'] = dirname(blastdb)
            blast_cmdline = NcbiblastpCommandline(
                db=basename(blastdb),
                outfmt=5,
                out=res_xml_fpath,
                num_threads=threads)
            out, err = blast_cmdline(stdin=str(rec.seq))
        else:
            res = blast_online(rec, res_xml_fpath)
            if res == 1:
                return 1

    search_for_best_hit(rec, res_xml_fpath, short_fpath)


blasted_singletones_dir = 'blasted_singletones'
new_proteomes_dir = 'new_proteomes'
new_annotations_dir = 'new_annotations'


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