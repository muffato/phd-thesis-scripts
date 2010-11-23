#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

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

