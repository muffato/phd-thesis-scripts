#! /users/ldog/muffato/python -OO

__doc__ = """
Lit sur l'entree standard des familles d'objets et les regroupe.
Il ne s'agit que d'une interface pour l'objet myCombinator.
"""

import sys
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs(["familiesFile"], [("showStats", bool, False)], __doc__)

# MAIN #

# On scanne toutes les lignes
f = utils.myTools.myOpenFile(noms_fichiers["familiesFile"], 'r')
comb = utils.myTools.myCombinator([])
for l in f:
	comb.addLink(l.split())
f.close()

# On affiche le resultat
for x in comb:
	print " ".join(x)
	if options["showStats"]:
		print >> sys.stderr, len(x)

