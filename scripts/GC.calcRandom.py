#! /users/ldog/muffato/python -OO

__doc__ = """
	Renvoie la difference moyenne entre 2 GC consecutifs
"""

import sys
import random
import utils.myTools

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["GCPercent"], [], __doc__)

# Chargement du fichier avec les taux de GC
f = utils.myTools.myOpenFile(noms_fichiers["GCPercent"], "r")
lst = []
dicGC = {}
for (j,ligne) in enumerate(f):
	gc = float(ligne.split()[8])
	lst.append(gc)
f.close()

for i in xrange(1000):
	random.shuffle(lst)
	random.shuffle(lst)
	random.shuffle(lst)
	count = 0.
	for j in xrange(len(lst)-1):
		count += abs(lst[j]-lst[j+1])
	print count/(len(lst)-1.)

