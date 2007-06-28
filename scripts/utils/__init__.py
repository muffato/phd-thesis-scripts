try:
	import utils.psyco
	utils.psyco.full()
	from utils.psyco.classes import *
	#import sys
	#print >> sys.stderr, "psyco fully loaded !"
except ImportError:
	import sys
	print >> sys.stderr, "Unable to load psyco !"

