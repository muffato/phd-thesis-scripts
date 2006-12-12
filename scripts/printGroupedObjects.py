#! /users/ldog/muffato/python -OO

__doc__ = """
Lit sur l'entree standard des familles d'objets et les regroupe.
Il ne s'agit que d'une interface pour l'objet myCombinator.
"""

import sys
import os
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs([], [("showStats", bool, False)], __doc__)

# MAIN #

comb = utils.myTools.myCombinator([])

# On scanne toutes les lignes
for l in sys.stdin:
	cc = l.split()
	comb.addLink(cc)

# On affiche le resultat
for x in comb.getGrp():
	print " ".join(x)
	if options["showStats"]:
		print >> sys.stderr, len(x)

