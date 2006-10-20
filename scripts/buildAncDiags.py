#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales de genes entre deux especes
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("ancetre",str,""), ("keepGaps",bool,True), ("fusionThreshold",int,-1), ("minimalLength",int,1), \
	("ancGenesFile",str,"/users/ldog/muffato/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers[1])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers[0], listEspeces)

genesAncRoot = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.root, False)
genesAncNode = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % options["ancetre"], False)

# 2. On genere toutes les diagonales entre des paires d'especes qui passent par le noeud

combin = utils.myTools.myCombinator([])

def combinDiag(c1, c2, nbTot, d):
	global combin, options
	if len(d) >= options["minimalLength"]:
		combin.addLink(d)

groupes = []
fils = set([])
for (e,_) in phylTree.items[options["ancetre"]]:
	tmp = phylTree.getSpecies(e)
	groupes.append(tmp)
	fils.update(tmp)
outgroup = set(listEspeces).difference(fils)

genomes = {}
#for e in fils:
for e in listEspeces:
	#genomes[e] = utils.myMaths.translateGenome(geneBank.dicEspeces[e], e, genesAncNode)
	genomes[e] = utils.myMaths.translateGenome(geneBank.dicEspeces[e], e, genesAncRoot)

#print groupes
for (i,j) in utils.myTools.myMatrixIterator(len(groupes), len(groupes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	for e1 in groupes[i]:
		for e2 in groupes[j]:
			print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
			utils.myMaths.iterateDiags(genomes[e1], genomes[e2], options["keepGaps"], options["fusionThreshold"], combinDiag)
			print >> sys.stderr, "OK"

#genomes = {}
#for e in listEspeces:
#	genomes[e] = utils.myMaths.translateGenome(geneBank.dicEspeces[e], e, genesAncRoot)

for g in groupes:
	for e1 in g:
		for e2 in outgroup:
			print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
			utils.myMaths.iterateDiags(genomes[e1], genomes[e2], options["keepGaps"], options["fusionThreshold"], combinDiag)
			print >> sys.stderr, "OK"

for g in combin.getGrp():
	print " ".join([str(x) for x in g])


