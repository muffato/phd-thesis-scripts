#! /usr/bin/python2.4

# Interface pour la classe myCombinator qui permet de regrouper des objets

import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myTools

(noms_fichiers, options) = myTools.checkArgs([], [("showStats", bool, False)], "Lit sur l'entree standard des familles d'objets et les regroupe")

# MAIN #

comb = myTools.myCombinator([])

# On scanne toutes les lignes
for l in sys.stdin:
	cc = l.split()
	comb.addLink(cc)

# On affiche le resultat
for x in comb.getGrp():
	print " ".join(x)
	if options["showStats"]:
		print >> sys.stderr, len(x)

