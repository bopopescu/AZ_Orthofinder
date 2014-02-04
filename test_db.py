#!/usr/bin/env python
from os.path import join, realpath

import traceback
from src.mysql import connector
from src.mysql.connector import errorcode

import sys
from site import addsitedir
addsitedir(join(realpath(__file__), 'src', 'mysql'))

try:
    cnx = connector.connect(
        user='orthomcl',
        password='1234',
        host='127.0.0.1',
        port='3307',
        database='orthomcl',
        buffered=True)
except:
    print 'Exception', traceback.print_exc(limit=0)
else:
    print 'Connected'