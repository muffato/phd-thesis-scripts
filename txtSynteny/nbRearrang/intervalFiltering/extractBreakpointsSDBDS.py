#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	A partir de la liste des etats des intervalles, cherche les points de rearrangements certains, de la forme Synt-Dup*-Rearr-Dup*-Synt
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( [("intervalStatus",file)], [], __doc__)

genome = {}
dicGenome = {}
nextChrom = 1

dicOldName = {}
status = {}
f = utils.myFile.openFile(arguments["intervalStatus"], "r")
for l in f:
	if l.startswith("-"):
		continue
	t = l.split()

	assert t[2] not in dicGenome
	if t[1] not in dicGenome:
		dicGenome[t[1]] = nextChrom
		genome[nextChrom] = [t[1]]
		nextChrom += 1
	c = dicGenome[t[1]]
	genome[c].append(t[2])
	dicGenome[t[2]] = c

	status[tuple(t[1:3])] = t[5:]
	dicOldName[t[1]] = t[3]
	dicOldName[t[2]] = t[4]
f.close()


def isRearr(x):
	return (status[x][1] in ["DUPLICATES_DIFF_ORIENT", "ONE_END", "NO_END", "DIFF_ORIENT"]) #and (status[x][0] == "WITHOUT_GENE_GAIN_INSIDE")

def isDup(x):
	return (status[x][1] == "DUPLICATES_SAME_ORIENT")

def isSynt(x):
	return (status[x][1] == "SAME_ORIENT")

for chrom in genome:

	ls = list(utils.myTools.myIterator.slidingTuple(genome[chrom]))

	for i in xrange(len(ls)):
		
		if not isRearr(ls[i]):
			continue

		j = i + 1
		while (j < len(ls)) and isDup(ls[j]):
			j += 1
		if not ((j < len(ls)) and isSynt(ls[j])):
			continue

		j = i - 1
		while (j >= 0) and isDup(ls[j]):
			j -= 1
		if not ((j >= 0) and isSynt(ls[j])):
			continue

		print utils.myFile.myTSV.printLine([ls[i][0], ls[i][1], dicOldName[ls[i][0]], dicOldName[ls[i][1]],  " ".join(status[ls[i]])])

