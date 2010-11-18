#!/usr/bin/env python2

__doc__ = """
	Compare les genes successifs dans le genome moderne et renvoie les scores de simultaneite
	Renvoie les memes scores pour les memes genes pris dans des ordres aleatoires
"""

import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree



# Score de ressemblance des profils des deux genes
def scoreSameEvents(v1, v2, lstindices):
	s = 0.
	n = 0
	for i in lstindices:
		s += v1[i] & v2[i]
		n += v1[i] | v2[i]
	return s/n if n != 0 else -1

@utils.myTools.memoize
def age(gene):
	for esp in reversed(ancGenesL):
		if gene in ancGenes[esp].dicGenes:
			return esp
	return arguments["modernSpecies"]

@utils.myTools.memoize
def nameInAnc(gene, esp):
	(c,i) = ancGenes[esp].dicGenes[gene]
	return ancGenes[esp].lstGenes[c][i].names[0]

def calcScore(lstGenes):
	s = []
	for i1 in xrange(len(lstGenes)):
		g1 = lstGenes[i1]
		a1 = age(g1)
		best = None
		for i2 in xrange(i1+1, len(lstGenes)):
			g2 = lstGenes[i2]
			a2 = age(g2)
			newcommon = a1 if phylTree.isChildOf(a1, a2) else a2
			if (best is None) or (phylTree.isChildOf(best, newcommon) and (best != newcommon)):
				indices = [phylTree.indNames[e] for e in phylTree.allNames if phylTree.isChildOf(e, newcommon) and (e != newcommon)]
				s1 = byname[nameInAnc(g1, newcommon)]
				s2 = byname[nameInAnc(g2, newcommon)]
				s.append( scoreSameEvents(s1, s2, indices) )
				best = newcommon
				if best == a1:
					break

	return s

	for (g1,g2) in utils.myTools.myIterator.slidingTuple(lstGenes):
		a1 = age(g1)
		a2 = age(g2)
		if phylTree.isChildOf(a1, a2):
			s1 = byname[nameInAnc(g1, a1)]
			s2 = byname[nameInAnc(g2, a1)]
			indices = [phylTree.indNames[e] for e in phylTree.allNames if phylTree.isChildOf(e, a1) and (e != a1)]
		else:
			s1 = byname[nameInAnc(g1, a2)]
			s2 = byname[nameInAnc(g2, a2)]
			indices = [phylTree.indNames[e] for e in phylTree.allNames if phylTree.isChildOf(e, a2) and (e != a2)]
		s.append( scoreSameEvents(s1, s2, indices) )
	return s

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("modernSpecies",str), ("ancestor",str), ("gene:treesSignatures",file)], [("modernGenomes",str,""), ("ancGenesFile",str,""), ("nbControl",int,1)], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
esp = arguments["modernSpecies"]
modernGenome = utils.myGenomes.Genome(arguments["modernGenomes"] % phylTree.fileName[esp])
ancGenes = phylTree.newCommonNamesMapperInstance()
ancGenes[arguments["modernSpecies"]] = modernGenome
ancGenesL = []
while esp != arguments["ancestor"]:
	(esp,_) = phylTree.parent[esp]
	ancGenes[esp] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[esp])
	ancGenesL.append(esp)

# Charge les signatures
class defaultdict2(dict):
	def __init__(self, val):
		self.defaultvalue = val

	def __getitem__(self, key):
		return dict.__getitem__(self, key) if key in self else self.defaultvalue

print >> sys.stderr, "Chargement des signatures ...",
allsign = {}
byname = defaultdict2([0] * len(phylTree.indNames))
f = utils.myFile.openFile(arguments["gene:treesSignatures"], "r")
for l in f:
	(name,_,val) = l.partition(' ')
	s = tuple(eval(val.replace('None', '0')))
	# Rassemble les signatures identiques -> sauve de la memoire
	byname[name] = allsign.setdefault(s, s)
f.close()
print >> sys.stderr, "OK"

genes = [g.names[0] for g in modernGenome]
for i in xrange(arguments["nbControl"]+1):
	l = calcScore(genes)
	for x in l:
		print i, x
	print >> sys.stderr, i, utils.myMaths.myStats.txtSummary(l)
	random.shuffle(genes)

