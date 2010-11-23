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

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("onlyOrthos",bool,False), ("IN.genesFile",str,""), ("IN.ancGenesFile",str,""), ("IN.diagsFile",str,"")], \
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
		diags[e].append([((c,i),l[i].strand) for i in xrange(len(l))])

for a in phylTree.listAncestr:
	# Les genes ancestraux
	genes[a] = utils.myGenomes.Genome(arguments["IN.ancGenesFile"] % phylTree.fileName[a])
	# Les diagonales
	diags[a] = []
	# On en profite pour lister les diagonales et les genes seuls
	notseen = set(xrange(len(genes[a].lstGenes[None])))
	f = utils.myFile.openFile(arguments["IN.diagsFile"] % phylTree.fileName[a], "r")
	for l in f:
		t = l.split("\t")
		d = [int(x) for x in t[2].split()]
		s = [int(x) for x in t[3].split()]
		s = [2*int(x>=0)-1 for x in s]
		diags[a].append(zip([(None,i) for i in d], s))
		notseen.difference_update(d)
	f.close()
	assert len(notseen) == 0
	#print >> sys.stderr, len(notseen)
	#diags[a].extend( [((None,g),0)] for g in notseen)

# Creation des dictionnaires genes -> diags
for (esp,lst) in diags.iteritems():
	dic = {}
	for (i,d) in enumerate(lst):
		for (j,(g,s)) in enumerate(d):
			dic[g] = (i,j,s)
	dicDiags[esp] = dic


# Les tables de conversion de diagonales entre ancetres successifs
val = {}
def do(node):

	# Les branches descendantes
	for (e,a) in phylTree.items.get(node, []):

		print >> sys.stderr, "Stating branch %s -> %s ..." % (node,e),

		# Les correspondances de genes entre les deux noeuds
		corresp = {}
		for (c,i) in dicDiags[e]:
			l = genes[node].getPosition(genes[e].lstGenes[c][i].names)
			assert len(l) in [0,1]
			if len(l) > 0:
				corresp[(c,i)] = dicDiags[node][tuple(l.pop())]
		print >> sys.stderr, 'e:', len(dicDiags[e]), 'node:', len(dicDiags[node]), 'corresp:', len(corresp),
		nbOK = 0
		nbNO = 0
		nbXX = 0
		for d in diags[e]:
			if arguments["onlyOrthos"]:
				d = [(g,s) for (g,s) in d if g in corresp]
			for ((g1,s1),(g2,s2)) in utils.myTools.myIterator.slidingTuple(d):
				if (g1 not in corresp) or (g2 not in corresp):
					continue
				(c1,i1,t1) = corresp[g1]
				(c2,i2,t2) = corresp[g2]
				j1 = i2 - s2*t2
				j2 = i1 + s1*t1
				print "XXX:", (s1,s2), corresp[g1], corresp[g2], (j1,j2), i2==j2, i1==j1, s1*s2 == t1*t2, len(diags[node][c1]), len(diags[node][c2]),
				if c1 == c2:
					if (i2 == j2) and (s1*s2 == t1*t2):
						nbOK += 1
						print "OK1"
					else:
						nbNO += 1
						print "NO1"
				elif (j2 >= 0) and (j2 < len(diags[node][c1])):
					print "NO2"
					nbNO += 1
				else:
					j1 = i2 - s2*t2
					if (j1 >= 0) and (j1 < len(diags[node][c2])):
						nbNO += 1
						print "NO3"
					else:
						nbXX += 1
						print "??"
		print node, e, nbOK, nbNO, nbXX, nbOK+nbNO+nbXX, (100.*nbOK)/(nbOK+nbNO), len(diags[node]), a
		print >> sys.stderr, "OK"
		val[(node,e)] = float(nbNO)/(nbOK+nbNO)
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


