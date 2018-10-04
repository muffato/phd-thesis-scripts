#!/usr/bin/env python2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
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
arguments = utils.myTools.checkArgs( \
	[("list",file), ("phylTree.conf",file), ("species",str)], \
	[("IN.ancGenesFile",str,""), ("IN.diagsFile",str,"")], \
	__doc__ \
)

# Chargement des tous les fichiers
###################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

ancGenes = {}
diags = {}


def areDup(names1, names2):
	anc = arguments["species"]
	while anc in phylTree.parent:
		(new,_) = phylTree.parent[anc]
		if new not in ancGenes:
			ancGenes[new] = utils.myGenomes.Genome(arguments["IN.ancGenesFile"] % new)

		g1 = ancGenes[new].getPosition(names1)
		g2 = ancGenes[new].getPosition(names2)
		if len(g1) == 0:
			return False
		if len(g2) == 0:
			return False
		assert len(g1) == 1
		assert len(g2) == 1
		g1 = g1.pop()[1]
		g2 = g2.pop()[1]

		if g1 == g2:
			return True

		names1 = ancGenes[new].lstGenes[None][g1].names
		names2 = ancGenes[new].lstGenes[None][g2].names
		anc = new
	
	return False



f = utils.myFile.openFile(arguments["list"], "r")
for l in f:
	l = l[:-1]
	t = l.split("\t")
	t = t[1].split("/")
	print "%s\t%d" % (l, areDup([t[0]], [t[1]]))
f.close()

