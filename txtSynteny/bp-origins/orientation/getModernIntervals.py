#!/usr/bin/env python2

__doc__ = """
	Prend deux listes d'intervalles et extrait les intervalles conserves 
"""

import sys
import math
import itertools
import collections

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( [("phylTree",file)], [("genes",str,""), ("ancGenes",str,"")], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])
phylTree.loadAllSpeciesSince(None, arguments["genes"])
anc = None

avgSize = {}
for e in phylTree.dicGenomes:
	lst = []
	genome = phylTree.dicGenomes[e]
	for c in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		for (g1,g2) in utils.myTools.myIterator.slidingTuple(genome.lstGenes[c]):
			lst.append(g2.beginning-g1.end-1)
	avgSize[e] = utils.myMaths.myStats.mean(lst)

for l in sys.stdin:
	t = l[:-1].split("\t")
	if t[0] != anc:
		ancGenes = utils.myGenomes.Genome(arguments["ancGenes"] % t[0])
		anc = t[0]
	g1 = int(t[1])
	g2 = int(t[2])
	for e in phylTree.species[anc]:
		genome = phylTree.dicGenomes[e]
		lg1 = genome.getPosition(ancGenes.lstGenes[None][g1].names)
		lg2 = genome.getPosition(ancGenes.lstGenes[None][g2].names)
		if (len(lg1) != 1) or (len(lg2) != 1):
			continue
		(c1,i1) = lg1.pop()
		(c2,i2) = lg2.pop()
		if c1 != c2:
			continue
		if abs(i1-i2) != 1:
			continue
		s1 = genome.lstGenes[c1][i1].strand
		s2 = genome.lstGenes[c2][i2].strand
		if i1 < i2:
			if (int(t[5]) != s1) or (int(t[6]) != s2):
				continue
			size = genome.lstGenes[c1][i2].beginning-genome.lstGenes[c1][i1].end-1
		else:
			if (int(t[5]) != -s1) or (int(t[6]) != -s2):
				continue
			size = genome.lstGenes[c1][i1].beginning-genome.lstGenes[c1][i2].end-1
		if size > 0:
			print utils.myFile.myTSV.printLine([ anc, e, t[7], math.log(size/avgSize[e], 2) ])

