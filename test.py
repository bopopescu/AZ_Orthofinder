from os import system, listdir, getcwd
import subprocess
import sys
prog_name = 'scenario_1.py'

for cmdline in \
    ['test_input/prots --prot-id-field 1',
     'test_input/annotations',
     'test_input/ids.txt -o test_input/test_ids',
     '--species-name test_input/species.txt -o test_input/test_ids']:

    res = subprocess.call('python ' + prog_name + ' ' + cmdline, shell=True)
    if res != 0:
        print >> sys.stderr, 'Command ' + prog_name + ' ' + cmdline + ' returned ' + str(res)
        exit(res)

    print
    print
    print '-' * 50