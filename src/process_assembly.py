from Bio import SeqIO
from Bio.Alphabet import IUPAC

import config
import logging
log = logging.getLogger(config.log_fname)


def filter_assembly(assembly, out, minimal_seq=1, skip=(), skip_after=None):
    records = []

    for i, rec in enumerate(SeqIO.parse(assembly, 'fasta')):
        if skip_after is not None and i > 10:
            continue
        #if len(rec.seq) <= minimal_seq:
        #    continue
        if i in skip:
            continue

        records.append(rec)

    if records:
        SeqIO.write(records, open(out, 'w'), 'fasta')
        return 0
    else:
        log.error('No sequences in %s are at least %d long.' % (assembly, minimal_seq))
        return 1