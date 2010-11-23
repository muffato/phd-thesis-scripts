#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import collections
import itertools

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
	diags[e] = []
	for (c,l) in genes[e].lstGenes.iteritems():
		diags[e].append([(c,i) for i in xrange(len(l))])

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
		diags[a].append([(None,i) for i in d])
	diags[a].extend( [(None,g)] for g in notseen)

# Creation des dictionnaires genes -> diags
for (esp,lst) in diags.iteritems():
	dic = {}
	for (i,d) in enumerate(lst):
		for (j,g) in enumerate(d):
			dic[g] = (i,j)
	dicDiags[esp] = dic


# Les tables de conversion de diagonales entre ancetres successifs
val = {}
def do(node):

	# Les branches descendantes
	for (e,l) in phylTree.items.get(node, []):

		print >> sys.stderr, "Stating branch %s -> %s ..." % (node,e),

		# Les correspondances de genes entre les deux noeuds
		corresp = {}
		for g in dicDiags[node]:
			corresp[g] = [dicDiags[e][tuple(x)] for x in genes[e].getPosition(genes[node].lstGenes[g[0]][g[1]].names)]

		nbTot = 0
		nbOK = 0
		for d in diags[node]:
			for (g1,g2) in utils.myTools.myIterator.slidingTuple(d):
				nbTot += 1
				for ((c1,i1),(c2,i2)) in itertools.product(corresp[g1], corresp[g2]):
					if (c1 == c2) and (abs(i1-i2) == 1):
						nbOK += 1
						break
		print node, e, nbOK, nbTot, (100.*nbOK)/nbTot, len(diags[node]), (100.*nbOK)/(nbTot+len(diags[node])), l
		print >> sys.stderr, "OK"
		val[(node,e)] = float(nbTot-nbOK)/nbTot
		do(e)

do(phylTree.root)

# Parcourt recursivement l'arbre et l'ecrit au format avec des parentheses, avec les longueurs de branche medianes
def convertToFlatFile(anc):

	a = phylTree.fileName[anc]
	if anc in phylTree.listSpecies:
		# On est arrive sur une feuille
		return a
	else:
		# On est sur un noeud, on construit la liste des distances
		l = []
		for (e,_) in phylTree.items[anc]:
			l.append(val[(anc,e)])
		# Constuction de la chaine finale
		return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for ((e,age),l) in zip(phylTree.items[anc],l) ]) + ")%s|%d" % (a,phylTree.ages[anc])


print convertToFlatFile(phylTree.root), ";"


