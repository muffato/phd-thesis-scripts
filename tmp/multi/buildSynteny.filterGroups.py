#!/usr/bin/env python2

__doc__ = """
	Separe les diagonales pairwise en 3 groupes en fonction de leur utilisation dans les blocs integres:
		1) celles utilisees a tous les ancetres intermediaires (= parfaites)
		2) celles non utilisees partout a cause de celles de 1) (= bloquees)
		3) le reste = celles non utilisees partout, mais a cause de 2) (= potentielles)
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import myDiags

# Arguments
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("target",str), ("diagGroups",file), ("ancGenomes",str)], [], __doc__ )


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les ancetres cibles
targets = set(myDiags.getTargets(phylTree, arguments["target"])[1])

# Chargement des blocs integres
ancGenomes = {}
for anc in targets:
	ancGenomes[anc] = utils.myGenomes.Genome(arguments["ancGenomes"] % phylTree.fileName[anc])

# Est-ce que la paire est vue dans l'ancetre
def isAncOK(anc, (g1,s1), (g2,s2)):
	(c1,i1) = ancGenomes[anc].dicGenes[g1]
	(c2,i2) = ancGenomes[anc].dicGenes[g2]
	if c1 != c2:
		return False
	if i2 != (i1 + s1*ancGenomes[anc].lstGenes[c1][i1].strand):
		return False
	if s1*ancGenomes[anc].lstGenes[c1][i1].strand != s2*ancGenomes[anc].lstGenes[c2][i2].strand:
		return False
	return True


# Remplissage des graphes et calcul du poids
f = utils.myFile.openFile(arguments["diagGroups"], "r")
for l in f:
	lst = eval(l)
	nbTot = 0
	nbOK = 0
	for x in lst:
		if x[0] in targets:
			nbTot += 1
			if isAncOK(x[0], (str(x[1]),x[3]), (str(x[2]),x[4])):
				nbOK += 1
	if nbTot > 0:
		print nbTot, nbOK, float(nbOK)/nbTot
	else:
		print >> sys.stderr, lst
f.close()

