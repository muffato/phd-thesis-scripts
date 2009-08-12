#! /users/ldog/muffato/python

__doc__ = """
	Lit un fichier et melange les lignes
"""

import random

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs([("fichier",file)], [("nbShuffle",int,10)], __doc__)

lst = utils.myFile.openFile(arguments["fichier"], 'r').readlines()
for _ in xrange(arguments["nbShuffle"]):
	random.shuffle(lst)

for x in lst:
	print x,

