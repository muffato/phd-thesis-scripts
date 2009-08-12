#! /users/ldog/muffato/python

__doc__ = """
	
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import os
import sys
import bisect
import operator
import itertools
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("diagsFile",file), ("ancGenesFile",file), ("modernGenesFile",file)], [], __doc__)

genome = utils.myGenomes.Genome(arguments["modernGenesFile"])
ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])

ancLabel = {}
for c in genome.lstGenes:
	ancLabel[c] = [None] * len(genome.lstGenes[c])

f = utils.myFile.myTSV.reader(arguments["diagsFile"])
for (n,t) in enumerate(f.csvobject):
	d = t[2].split()
	if len(d) == 1:
		continue
	for x in d:
		x = int(x)
		for (c,i) in genome.getPosition(ancGenes.lstGenes[None][x].names):
			ancLabel[c][i] = n


f.file.close()

for c in genome.lstChr:
	lastn = None
	laststart = None
	lastend = None

	for (g,n) in itertools.izip(genome.lstGenes[c], ancLabel[c]):
		#print c, g.names[0], g.beginning, g.end, n
		if n is None:
			continue
		if n != lastn:
			if laststart is not None:
				#if lastend != laststart:
					print c, laststart.beginning, lastend.end, lastn
			laststart = g
			lastn = n
		lastend = g
	



