
import sys

if "+psyco" in sys.argv:
	sys.argv.remove("+psyco")
	try:
		import utils.psyco
		utils.psyco.full()
		from utils.psyco.classes import __metaclass__
	except ImportError:
		print >> sys.stderr, "Unable to load psyco !"

