from ftplib import FTP
from genericpath import isdir
from os import mkdir, makedirs, remove
from os.path import join, basename
from Bio import SeqIO

species = ['Klebsiella pneumoniae', 'Escherichia coli']

proteomes_dir = '../../EcoliKpneumoniae/proteomes'
if not isdir(proteomes_dir):
    makedirs(proteomes_dir)

ftp = FTP('ftp.ncbi.nih.gov')
print ftp.login()
ftp.cwd('genomes/Bacteria')
for dirname in ftp.nlst():
    for sp in species:
        sp = sp.replace(' ', '_')
        if dirname.startswith(sp):
            print dirname
            for fpath in ftp.nlst(dirname):
                if fpath.endswith('.gbk'):
                    print fpath
                    dest = join(proteomes_dir, basename(fpath))
                    ftp.retrbinary('RETR ' + fpath, open(dest, 'wb').write)

                    rec = SeqIO.read(dest, 'gb')
                    print '       Definition: ' + rec.description
                    if 'plasmid' in rec.description:
                        remove(dest)
                    else:
                        print '       saved ' + dest
                    print

                    #rec = SeqIO.read(gb_fpath, genbank_ext)
                    #log.info('       Organism: ' + rec.annotations['organism'])
                    #log.info('       Definition: ' + rec.description)
                    #if 'plasmid' in rec.description:
                    #    remove(gb_fpath)
                    #else:
                    #    log.info('       saved %s' % gb_fpath)
                    #    print '   ' + fname

ftp.quit()