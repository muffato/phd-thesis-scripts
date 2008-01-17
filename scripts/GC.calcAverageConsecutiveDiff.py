#! /users/ldog/muffato/python -OO

__doc__ = """
	Renvoie la difference moyenne entre 2 GC consecutifs
"""

import sys
import random
import utils.myMaths
import utils.myTools

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["GCPercent"], [("GCcolumn",int,0), ("nbShuffle",int,3)], __doc__)

# Chargement du fichier avec les taux de GC
f = utils.myTools.myOpenFile(noms_fichiers["GCPercent"], "r")
lst = []
for (j,ligne) in enumerate(f):
	gc = float(ligne.split()[options["GCcolumn"]])
	lst.append(gc)
f.close()

for _ in xrange(options["nbShuffle"]):
	random.shuffle(lst)

count = []
for (gc1,gc2) in utils.myTools.slidingTuple(lst):
	count.append(abs(gc1-gc2))
print utils.myMaths.myStats(lst)

