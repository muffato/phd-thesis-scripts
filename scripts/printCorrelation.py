#! /users/ldog/muffato/python -OO

import sys
import utils.myTools
import utils.myMaths

(fichiers,options) = utils.myTools.checkArgs(["dataFile"], [("column1",int,0), ("column2",int,1)], "Renvoie la correlation entre les deux jeux de donnees")

lx1 = []
lx2 = []

f = utils.myTools.myOpenFile(fichiers["dataFile"], 'r')
for l in f:
	t = l.replace('\n', '').split("\t")
	try:
		d1 = float(t[options["column1"]])
		d2 = float(t[options["column2"]])
		lx1.append(d1)
		lx2.append(d2)
	except ValueError:
		pass

f.close()

"""print >> sys.stderr, "Loading ...",
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

print >> sys.stderr, len(lx1), "OK"""

print utils.myMaths.correlation(lx1, lx2)

