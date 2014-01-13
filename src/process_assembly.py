from Bio import SeqIO

import config
import logging
log = logging.getLogger(config.log_fname)


def filter_assembly(assembly, out):
    records = []

    for rec in SeqIO.parse(assembly, 'fasta'):
        if len(rec.seq) > 200000:
            records.append(rec)

    if records:
        SeqIO.write(records, out, 'fasta')
        return 0
    else:
        log.error('No seqences in %s are at least 20k long.' % assembly)
        return 1