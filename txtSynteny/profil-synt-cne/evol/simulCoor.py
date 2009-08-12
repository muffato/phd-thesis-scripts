#! /users/ldog/muffato/python

__doc__ = """
	Simule les points de cassure du genome en fonction des tailles des espaces
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
arguments = utils.myTools.checkArgs( [("refGenome",file)], [], __doc__)

genome = utils.myGenomes.Genome(arguments["refGenome"])

ages = [0,5,15,30,55,60,90,95,105,130,175,320,340,420,510,550,580]
nb = [19525,16844,15587,14728,14304,14220,13980,13689,11251,10487,8309,7585,5670,3907,241,231,9]


allInt = set()
for c in genome.lstGenes:
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(genome.lstGenes[c]):
		if g1.end < g2.beginning:
			allInt.add( ((g2.beginning-1)-(g1.end+1)+1, g1, g2) )
allIntL = sorted(allInt)
rand = utils.myMaths.randomValue()
#for i in xrange(20000):
while True:

	#if (i % 100) == 0:
	if len(allIntL) == nb[0]:
		count = collections.defaultdict(int)
		for (d,g1,g2) in allIntL:
			print ages[0], d
			#print i, d, g1.strand, g2.strand
			count[(g1.strand,g2.strand)] += 1
		s = float(sum(count.values()))/100
		#print ages[0], (count[(1,1)]+count[(-1,-1)])/s, count[(1,-1)]/s, count[(-1,1)]/s
		ages.pop(0)
		nb.pop(0)
	sizes = [x[0] for x in allIntL]
	rand.initBisect(sizes)
	pos = rand.getRandomBisectPos()
	#(l,g1,g2) = allIntL[pos]
	del allIntL[pos]

