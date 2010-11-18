#!/usr/bin/env python2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( \
	[("list",file), ("phylTree.conf",file), ("species",str)], \
	[("IN.ancGenesFile",str,""), ("IN.diagsFile",str,"")], \
	__doc__ \
)

# Chargement des tous les fichiers
###################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

diags = {}
ancGenes = {}
anc = arguments["species"]
while anc in phylTree.parent:
	(new,_) = phylTree.parent[anc]
	ancGenes[new] = utils.myGenomes.Genome(arguments["IN.ancGenesFile"] % new)

	dicDiags = {}
	for (i,l) in enumerate(utils.myFile.myTSV.readTabular(arguments["IN.diagsFile"] % new, [str,int,(str,4)])):
		#d = [int(x) for x in l[2].split()]
		#s = [int(x) for x in l[3].split()]
		for (j,(g,s)) in enumerate(itertools.izip(l[2].split(),l[3].split())):
			dicDiags[int(g)] = (i, j, int(s))
	diags[new] = dicDiags

	anc = new

def getGeneAge(names):
	age = 0
	for (anc,genes) in ancGenes.iteritems():
		g = genes.getPosition(names)
		if len(g) != 0:
			age = max(phylTree.ages[anc], age)
	return age
	

def getIntergenicAge(names1, names2, strand1, strand2):

	age = 0
	for (anc,dicDiags) in diags.iteritems():
		
		g1 = ancGenes[anc].getPosition(names1)
		g2 = ancGenes[anc].getPosition(names2)
		
		if (len(g1) > 0) and (len(g2) > 0):
			assert len(g1) == 1
			assert len(g2) == 1
			g1 = g1.pop()[1]
			g2 = g2.pop()[1]

			(c1,i1,s1) = dicDiags[g1]
			(c2,i2,s2) = dicDiags[g2]

			if (c1 == c2) and ((i2==i1+1 and s1==strand1 and s2==strand2) or (i2==i1-1 and s1==-strand1 and s2==-strand2)):
				age = max(age, phylTree.ages[anc])
	return age


f = utils.myFile.openFile(arguments["list"], "r")
for l in f:
	l = l[:-1]
	t = l.split("\t")
	if t[0] == "GENE":
		print l + "\t" + str(getGeneAge([t[1]]))
	else:
		s = t[5].split("/")
		t = t[1].split("/")
		print l + "\t" + str(getIntergenicAge([t[0]], [t[1]], int(s[0]), int(s[1])))
f.close()

