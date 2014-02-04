import sys
import traceback
import src.mysql.connector
from src.mysql.connector import errorcode

try:
    cnx = src.mysql.connector.connect(
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