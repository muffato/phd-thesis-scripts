#! /users/ldog/muffato/python -OO

import sys
import utils.myTools
import utils.myMaths

(fichiers,_) = utils.myTools.checkArgs(["fichier1", "fichier2"],[],"Renvoie la correlation entre les deux jeux de donnees")

print >> sys.stderr, "Loading ...",
f1 = utils.myTools.myOpenFile(fichiers["fichier1"], 'r')
f2 = utils.myTools.myOpenFile(fichiers["fichier2"], 'r')
lx1 = []
lx2 = []
for l1 in f1:
	l2 = f2.next()
	try:
		x1 = float(l1)
		x2 = float(l2)
		lx1.append(x1)
		lx2.append(x2)
	except ValueError:
		pass

print >> sys.stderr, len(lx1), "OK"

print utils.myMaths.correlation(lx1, lx2)

