#! /users/ldog/muffato/python

import sys
import utils.myTools
import utils.myMaths

arguments = utils.myTools.checkArgs([("dataFile",file)], [("column1",int,0), ("column2",int,1)], "Renvoie la correlation entre les deux jeux de donnees")

lx1 = []
lx2 = []

f = utils.myTools.myOpenFile(arguments["dataFile"], 'r')
for l in f:
	t = l.replace('\n', '').split("\t")
	try:
		d1 = float(t[arguments["column1"]])
		d2 = float(t[arguments["column2"]])
		lx1.append(d1)
		lx2.append(d2)
	except ValueError:
		pass

f.close()

print utils.myMaths.myStats.correlation(lx1, lx2)

