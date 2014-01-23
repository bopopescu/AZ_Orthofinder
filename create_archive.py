#!/usr/bin/env python

from genericpath import isfile, exists
from os import walk, remove, mkdir, symlink, unlink
from os.path import isdir, join, relpath, abspath
from shutil import copytree, ignore_patterns, copy2, rmtree, copy
from zipfile import ZipFile, ZIP_DEFLATED
from datetime import datetime
from time import strftime, strptime

about_html = 'about.html'

structure = [
    about_html,
    'scenario_1.py',
    'scenario_2.py',
    'clean_db.py',
    'src',
    '.cpan',
    'test_input',
    'config.txt',
    ]

dropbox_folder = '/Users/vladsaveliev/Dropbox/Public/az'

archive_dir = 'orthofinder'

ignored = '*.pyc', 'CVS', '.git', 'tmp', '.svn', '.DS_Store', '*.o'

def compress(source_dir, dest_file):
    relroot = abspath(join(source_dir, ".."))

    with ZipFile(dest_file, 'w', compression=ZIP_DEFLATED) as zf:
        for root, dirs, files in walk(source_dir):
            zf.write(root, relpath(root, relroot))

            for f in files:
                fname = join(root, f)
                if isfile(fname):
                    arcname = join(relpath(root, relroot), f)
                    #if '.cpan' not in fname and 'perl_modules' not in fname:
                    print '  ' + arcname  # join(relpath(root, directory), f)
                    zf.write(fname, arcname)

if __name__ == '__main__':
    if isdir(archive_dir):
        print 'Error: directory %s exists' % archive_dir
        exit(1)

    mkdir(archive_dir)
    print 'Created temporary directory %s' % archive_dir

    for obj in structure:
        if isfile(obj):
            copy2(obj, join(archive_dir, obj))
        elif isdir(obj):
            copytree(obj, join(archive_dir, obj), ignore=ignore_patterns(*ignored))
        else:
            print '  warning: %s does not exist.' % obj

    archive_basename = archive_dir + '_' + \
                       datetime.now().strftime('%Y_%m_%d')
    i = 1
    archive_file = archive_basename + '.zip'
    while exists(archive_file):
        i += 1
        archive_file = archive_basename + '__' + str(i) + '.zip'

    print '\nCompessing to %s' % archive_file
    compress(archive_dir, archive_file)

    print
    print 'Copying %s to %s' % (about_html, dropbox_folder)
    if exists(join(dropbox_folder, about_html)):
        remove(join(dropbox_folder, about_html))
    copy2(about_html, dropbox_folder)

    if exists(join(dropbox_folder, archive_file)):
        remove(join(dropbox_folder, archive_file))
    print 'Copying %s to %s' % (archive_file, dropbox_folder)
    copy2(archive_file, dropbox_folder)

    rmtree(archive_dir)

    copy_fpath = join(dropbox_folder, archive_dir + '.zip')
    if isfile(copy_fpath):
        remove(copy_fpath)
    copy(archive_file, copy_fpath)
    #try:
    #    unlink(symlink_fpath)
    #except:
    #    pass
    #symlink(archive_file, symlink_fpath)