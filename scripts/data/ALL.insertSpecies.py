#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les arbres de proteines et rajoute une nouvelle espece a partir d'orthologues
"""

import os
import sys
import time
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

# Arguments
############
(noms_fichiers, options) = utils.myTools.checkArgs( ["newPhylTree.conf", "proteinTree", "orthologues"], [("newSpecies",str,""), ("oldAncGenes",str,"")], __doc__ )

# Chargement des fichiers
##########################
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["newPhylTree.conf"])
orthologues = utils.myGenomes.Genome(noms_fichiers["orthologues"])
(data,info,roots) = utils.myProteinTree.loadTree(noms_fichiers["proteinTree"])



####################################################################
# On considere que les duplications 'dubious' ne sont pas valables #
####################################################################
def isDuplicatedNode(inf):
	return (inf['Duplication'] != 0) and ('dubious_duplication' not in inf)


######################################################################
# Construction du dictionnaire: nom_proteine -> noeud/arbre/nom_gene #
######################################################################
def store(node):
	if node in data:
		for (fils,_) in data[node]:
			store(fils)
	else:
		dicNames[info[node]['protein_name']] = (i,node,info[node]['gene_name'])


##############################################
# Construction du dictionnaire: fils -> pere #
##############################################
def mkPar(node):
	if node in data:
		for (fils,_) in data[node]:
			dicPar[fils] = node
			mkPar(fils)


# Initialisation
#################
maxNodeID = max(info)
(newNode,age0) = phylTree.parent[options["newSpecies"]]

dicNames = {}
dicTrees = {}
for (i,r) in enumerate(roots):
	dicTrees[i] = r
	store(r)

# On definit les endroits ou rajouter les familles
###################################################
for famille in orthologues:

	# Les genes des anciennes especes
	lstPar = [dicNames[g] for g in famille.names if g in dicNames]
	
	# On verifie que les liens ne menent qu'a un seul arbre
	ntrees = len(set([x[0] for x in lstPar]))
	if ntrees != 1:
		print >> sys.stderr, "> PB NB FAM", ntrees
		continue
	tree = dicTrees[lstPar[0][0]]

	# 1. Identifier l'ancetre commun

	# La table des parents
	dicPar = {}
	mkPar(tree)
	names = [x[1] for x in lstPar]

	# On remonte du gene de reference a la racine
	n0 = n = names.pop()
	refPath = {}
	while n != tree:
		n = dicPar[n]
		refPath[n] = (len(refPath),n)

	# On remonte les autres genes en conservant celui qui va le plus haut
	best = (-1,n0)
	while len(names) > 0:
		n = names.pop()
		while n not in refPath:
			n = dicPar[n]
		best = max(best, refPath[n])

	# -> ancestor = Le noeud de l'ancetre commun
	ancestor = best[1]

	# 2. Identifier le point d'insertion

	# On va remonter l'arbre jusqu'a ce qu'on depasse la cible
	path = []
	while phylTree.ages.get(info[ancestor]['taxon_name']) < age0:
		path.append(ancestor)
		if ancestor in dicPar:
			ancestor= dicPar[ancestor]
		else:
			print >> sys.stderr, "NEW ROOT", ancestor
			break
	else:
		# Puis on redescend pour eviter les duplications
		while isDuplicatedNode(info[ancestor]):
			ancestor = path.pop()
		print >> sys.stderr, "INSERTION AFTER", ancestor
	
	# 3. Insertion
	
	# Creer les feuilles de la nouvelle espece
	maxNodeID += 1
	esp = []
	for g in famille.names:
		if g in dicNames:
			maxNodeID += 1
			info[maxNodeID] = {'taxon_name':options["newSpecies"], 'Duplication':0, 'Bootstrap':-1, 'gene_name':g}
			esp.append(maxNodeID)
	# Fusions que l'arbre soit binaire
	while len(esp) >= 2:
		e1 = esp.pop()
		e2 = esp.pop()
		maxNodeID += 1
		info[maxNodeID] = {'taxon_name':options["newSpecies"], 'Duplication':1, 'Bootstrap':-1}
		data[maxNodeID] = [(e1,0), (e2,0)]
		esp.append(maxNodeID)
	espNode = maxNodeID

	# Il faut creer le noeud juste au dessus de ancestor
	maxNodeID += 1
	info[maxNodeID] = {'taxon_name':newNode, 'Duplication':0, 'Bootstrap':-1}
	data[maxNodeID] = [(ancestor,0),(espNode,0)]

	if ancestor in dicPar:
		p = dicPar[ancestor]
		newData = []
		for x in data[p]:
			if x[0] == ancestor:
				newData.append( (maxNodeID,x[1]) )
			else:
				newData.append( x )
		data[p] = newData
		# Quitte a renommer les noeuds de duplication
		while isDuplicatedNode(info[p]):
			if phylTree.ages.get(info[p]['taxon_name']) > age0:
				break
			info[p]['taxon_name'] = newNode
			p = dicPar[p]
	else:
		dicTrees[lstPar[0][0]] = maxNodeID


print >> sys.stderr, "Mise en forme des arbres ...",
nb = 0
for r in dicTrees.itervalues():
	nb += 1
	utils.myProteinTree.printTree(sys.stdout, data, info, r)
print >> sys.stderr, "%d arbres OK" % nb


