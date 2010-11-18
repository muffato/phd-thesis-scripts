#!/usr/bin/env python2

__doc__ = """
	Lit des listes d'id de genes et renvoie les noms a la place
"""


import sys

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("genesAncestraux",file), ("input",file)], [], __doc__)

genesAnc = utils.myGenomes.Genome(arguments["genesAncestraux"])

f = utils.myFile.openFile(arguments["input"], "r")
for l in f:
	t = l.split()
	for gs in t:
		g = genesAnc.lstGenes[None][int(gs)]
		print g.names[0],
	print
f.close()

