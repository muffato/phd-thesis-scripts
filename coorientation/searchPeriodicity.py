#! /users/ldog/muffato/python -OO

__doc__ = """
Lit un fichier de distances et en selectionne un certain nombre, selon un modele neutre de cassures proportionnelles aux distances
"""

import sys
import utils.myMaths
import utils.myTools


(noms_fichiers,options) = utils.myTools.checkArgs(["distFile"], [("binSize",int,100)], __doc__)

print >> sys.stderr, "Chargement des distances ...",
f = utils.myTools.myOpenFile(noms_fichiers["distFile"], "r")
allDist = []
allCoorientation = []
for ligne in f:
	t = ligne.split()
	allDist.append(int(t[1]))
	allCoorientation.append(intern(t[0]))
f.close()
print >> sys.stderr, len(allDist), "OK"

print >> sys.stderr, "Mise en paquets de taille %d ..." % options["binSize"],
bins = [0] * (1 + max(allDist) / options["binSize"])
for d in allDist:
	bins[d/options["binSize"]] += 1
print >> sys.stderr, sum(bins), "OK"

def autocorrel(n):
	s = 0
	for i in xrange(len(bins)-n):
		s += bins[i]*bins[i+n]
	return s / (len(bins)-n)

allCorrel = [autocorrel(i) for i in xrange(1,len(bins)-1)]

tmp = [allCorrel.index(x) for x in sorted(allCorrel, reverse=True)[:100]]

print tmp
print max(tmp) * options["binSize"]

