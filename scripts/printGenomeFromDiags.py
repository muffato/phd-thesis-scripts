#! /users/ldog/muffato/python

__doc__ = """
Convertit un genome (suite de diagonales) en genome (suite de genes)
"""

import sys
import utils.myTools
import utils.myGenomes

itStrictUpperMatrix = utils.myTools.myIterator.tupleOnStrictUpperList

arguments = utils.myTools.checkArgs( [("communitiesFile",file), ("ancGenesFile",file)], [], __doc__)

genes = [g for g in utils.myGenomes.Genome(arguments["ancGenesFile"])]

lastC = None
f = utils.myTools.myOpenFile(arguments["communitiesFile"], 'r')
for l in f:
	t = l.replace('\n', '').split('\t')
	if t[0] != lastC:
		lastC = t[0]
		nbC = 0
	diag = [int(x) for x in t[-2].split()]
	strand = [int(x) for x in t[-1].split()]
	if len(t) == 4:
		if int(t[1]) < 0:
			diag.reverse()
			strand = [-x for x in strand.__reversed__()]
	for (i,g) in enumerate(diag):
		print utils.myFile.myTSV.printLine([lastC, nbC+i, nbC+i+1, strand[i], " ".join(genes[g].names)])
	nbC += len(diag)

