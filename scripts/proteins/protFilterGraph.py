#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit un graphe (aretes + noms) et ne selectionne que les noeuds desires
	  Produit un nouveau graphe (aretes + noms)
"""

import sys
import math
import utils.myTools



arguments = utils.myTools.checkArgs( [("grapheInit",file), ("lstNoms",file)], [], __doc__ )

# MAIN #

# On lit la liste des fusions
lstOK = set()
f = utils.myTools.myOpenFile(arguments["lstNoms"], 'r')
for l in f:
	lstOK.update(l.split())

# Le nouveau graphe
fg = utils.myTools.myOpenFile(arguments["grapheInit"], 'r')
nb = 0
for l in fg:
	
	t = l.split()

	if (t[0] in lstOK) and (t[1] in lstOK):
		print l,
		nb += 1
fg.close()

print >> sys.stderr, "%d aretes gardees" % nb

