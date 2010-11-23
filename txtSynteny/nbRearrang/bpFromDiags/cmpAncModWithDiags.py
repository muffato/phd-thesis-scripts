#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Parcourt les deux genomes et extrait les intervalles conserves, les points de cassure
	  en tenant compte des evenements de genes (gain/perte/duplication)
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
arguments = utils.myTools.checkArgs( [("ancGenome",file), ("modernGenome",file), ("diags",file)], [], __doc__)

ancGenome = utils.myGenomes.Genome(arguments["ancGenome"])
genome = utils.myGenomes.Genome(arguments["modernGenome"])

dicA = {}
dicM = {}
dicD = {}

f = utils.myFile.openFile(arguments["diags"], "r")
for (i,l) in enumerate(f):
	t = l[:-1].split("\t")
	dS = [int(x) for x in t[9].split()]
	dA = t[4].split()
	dM = t[7].split()
	for (x,s) in itertools.izip(dA, dS):
		dicA[x] = (i,s)
	for (x,s) in itertools.izip(dM, dS):
		dicM[x] = (i,s)
	dicD[i] = (dA, dM)
f.close()
print >> sys.stderr, dicD

def rewriteGenome(genome, dic):
	newGenome = {}
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		tmp = [(dic[gene.names[0]],gene.strand) for gene in genome.lstGenes[chrom] if gene.names[0] in dic]
		tmp = [(i,s*sa) for ((i,sa),s) in tmp]
		if len(tmp) > 0:
			newGenome[chrom] = [i for (i,l) in itertools.groupby(tmp)]
	return newGenome

newGA = rewriteGenome(ancGenome, dicA)
newGM = rewriteGenome(genome, dicM)

for cA in newGA:
	print >> sys.stderr, "ANC", cA, newGA[cA]
for cM in newGM:
	print >> sys.stderr, "MOD", cM, newGM[cM]

def getAdj(genome):
	adj = set()
	extr1 = set()
	extr2 = set()
	for (chrom,l) in genome.iteritems():
		adj.update(utils.myTools.myIterator.slidingTuple(l))
		extr1.add( l[0] )
		extr2.add( l[-1] )
		extr1.add( (l[-1][0],-l[-1][1]) )
		extr2.add( (l[0][0],-l[0][1]) )
	return (adj,extr1,extr2)

(adjA,extrA1,extrA2) = getAdj(newGA)
(adjM,extrM1,extrM2) = getAdj(newGM)

def check(adjA, adjM, extrM1, extrM2, txt):
	def span(l):
		return "%s >%d> %s" % (l[0], len(l), l[-1])
	for ((i,si),(j,sj)) in adjA:
		if (((i,si),(j,sj)) not in adjM) and (((j,-sj),(i,-si)) not in adjM):
			if ((i,si) not in extrM2) or ((j,sj) not in extrM1):
				#print txt, (i,j), dicD[i], dicD[j]
				print txt, "%d/%d" % (i,si), "+", "%d/%d" % (j,sj), ":", span(dicD[i][0]), "/", span(dicD[i][1]), "!=!", span(dicD[j][0]), "/", span(dicD[j][1])

check(adjA, adjM, extrM1, extrM2, "LOST")
check(adjM, adjA, extrA1, extrA2, "GAIN")

