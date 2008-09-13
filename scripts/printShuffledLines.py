#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit un fichier et melange les lignes
"""

import random
import utils.myTools


arguments = utils.myTools.checkArgs([("fichier",file)], [("nbShuffle",int,10)], __doc__)

lst = utils.myTools.myOpenFile(arguments["fichier"], 'r').readlines()
for _ in xrange(arguments["nbShuffle"]):
	random.shuffle(lst)

for x in lst:
	print x,

