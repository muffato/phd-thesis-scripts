#! /users/ldog/muffato/python -OO

__doc__ = """
Lit un fichier de distances et en selectionne un certain nombre, selon un modele neutre de cassures proportionnelles aux distances
"""

import sys
import utils.myMaths
import utils.myTools


(noms_fichiers,options) = utils.myTools.checkArgs(["distFile"], [("nbConserved",int,0)], __doc__)

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

print >> sys.stderr, "Simulation ...",
notChosen = set(xrange(len(allDist)))
randomPicker = utils.myMaths.randomValue(allDist)
nb = 0
while len(notChosen) > options["nbConserved"]:
	nb += 1
	notChosen.discard(randomPicker.getRandomPos())
print >> sys.stderr, nb, "OK"

for x in notChosen:
	print allCoorientation[x], allDist[x]
