#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les arbres de proteines et rajoute la nouvelle espece aves les orthologues
"""

import os
import sys
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( ["newPhylTree.conf", "proteinTree", "orthologues"], [("newNode",str,""), ("oldAncGenes",str,"")], __doc__ )

# Chargement des fichiers
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["newPhylTree.conf"])
orthologues = utils.myGenomes.Genome(noms_fichiers["orthologues"])
(data,info,roots) = utils.myProteinTree.loadTree(noms_fichiers["proteinTree"])

par = phylTree.parent[options["newNode"]][0]

# On construit un dictionnaire "proteines" -> "genes"
dicNames = {}
def store(node):
	pofficialName[node] = node
	if node in data:
		for (fils,length) in data[node]:
			store(fils)
	else:
		dicNames[info[node]['protein_name']] = (r,node,info[node]['gene_name'])
pofficialName = {}
for r in roots:
	store(r)


print >> sys.stderr, "go go go !"
# On definit les endroits ou rajouter les familles
allPhyl = {}
age0 = phylTree.ages[par]
for famille in orthologues:
	new = [g for g in famille.names if g not in dicNames]
	lstPar = [dicNames[g] for g in famille.names if g in dicNames]
	ntrees = len(set([x[0] for x in lstPar]))
	if ntrees > 1:
		print "> MULTI"
		continue
	elif ntrees == 0:
		print "> NO FAMILY"
		continue
	r = lstPar[0][0]
	tree = utils.myPhylTree.PhylogeneticTree( (data,r,{}) )
	tree.officialName = pofficialName
	print "> LAST COMMON ANCESTOR"
	print lstPar
	ii = tree.lastCommonAncestor([x[1] for x in lstPar])
	print ii
	print info[ii]
	# On va remonter l'arbre jusqu'a ce qu'on arrive avant le pere
	while phylTree.ages.get(info[ii]['taxon_name']) < age0:
		if ii in tree.parent:
			(ii,_) = tree.parent.get(ii)
		else:
			print "NEW ROOT", ii
			break
		print info[ii]
	else:
		print "INSERTION AFTER", ii
	print



sys.exit(0)

print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
geneFamilies = utils.myTools.defaultdict(list)
for r in roots:
	nb += 1
	extractGeneFamilies(r, "FAM%d" % nb, None, None)
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, "%d arbres OK" % nb


