from site import addsitedir
from os.path import join, realpath
addsitedir(join(realpath(__file__), '..', '..'))