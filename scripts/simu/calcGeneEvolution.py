#! /users/ldog/muffato/python

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("genesFile",str), ("ancGenesFile",str)], [], __doc__)

# Chargement des tous les fichiers
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
genes = {}
for e in phylTree.listSpecies:
	genes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
	genes[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])

def transformName(esp, (c,i)):
	if esp in phylTree.items:
		return i
	else:
		return str(c) + "|" + str(i)

def do(node):

	for (e,_) in phylTree.items.get(node, []):
		res = []
		seen = set([transformName(e,(c,i)) for (c,l) in genes[e].lstGenes.iteritems() for i in xrange(len(l))])
		for g in genes[node].lstGenes[None]:
			lnewg = [transformName(e,x) for x in genes[e].getPosition(g.names)]
			seen.difference_update(lnewg)
			res.append(lnewg)
		print (res,seen)
		do(e)

print (phylTree.root,len(genes[phylTree.root].lstGenes[None]))

do(phylTree.root)

