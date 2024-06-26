#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import collections

import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.genesFile",str,""), ("IN.ancGenesFile",str,""), ("IN.diagsFile",str,"")], \
	__doc__ \
)

# Chargement des tous les fichiers
###################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

genes = {}
diags = {}
dicDiags = {}

for e in phylTree.listSpecies:
	# Les genes des especes modernes
	genes[e] = utils.myGenomes.Genome(arguments["IN.genesFile"] % phylTree.fileName[e])
	diags[e] = [[g] for g in xrange(len(list(genes[e])))]

for a in phylTree.listAncestr:
	# Les genes ancestraux
	genes[a] = utils.myGenomes.Genome(arguments["IN.ancGenesFile"] % phylTree.fileName[a])
	# Les diagonales
	tmp = utils.myDiags.loadDiagsFile(arguments["IN.diagsFile"] % phylTree.fileName[a], [a], phylTree.officialName)[a]
	# On en profite pour lister les diagonales et les genes seuls
	notseen = set(xrange(len(genes[a].lstGenes[None])))
	diags[a] = []
	for (d,_,_,_,_) in tmp:
		notseen.difference_update(d)
		diags[a].append(d)
	diags[a].extend( [g] for g in notseen)

# Creation des dictionnaires genes -> diags
for (esp,lst) in diags.iteritems():
	dic = {}
	for (i,d) in enumerate(lst):
		for (j,g) in enumerate(d):
			dic[g] = (i,j)
	dicDiags[esp] = dic
	print (esp,lst)


# Les tables de conversion de diagonales entre ancetres successifs
def do(node):

	# Les branches descendantes
	for (e,_) in phylTree.items.get(node, []):

		print >> sys.stderr, "Stating branch %s -> %s ..." % (node,e),

		# Les correspondances de genes entre les deux noeuds
		corresp = [[j for (_,j) in genes[node].getPosition(g.names)] for g in genes[e]]

		res = []
		for d in diags[e]:
			count = collections.defaultdict(int)
			for g in d:
				for gg in corresp[g]:
					if gg in dicDiags[node]:
						count[dicDiags[node][gg][0]] += 1
			res.append(count.items())
		print res
		print >> sys.stderr, "OK"
		do(e)

do(phylTree.root)

