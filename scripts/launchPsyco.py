#! /usr/bin/python2.5 -OO

import sys

try:
	import utils.psyco
	utils.psyco.full()
	from utils.psyco.classes import *
except ImportError:
	pass

sys.argv = sys.argv[1:]
execfile(sys.argv[0])

