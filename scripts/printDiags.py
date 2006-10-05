#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales de genes chaque couple d'especes.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


########
# MAIN #
########


# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["geneList.conf", "orthologuesList"], \
	[("speciesUsed",str,""), ("keepGaps",bool,True), ("fusionThreshold",int,-1)], \
	__doc__ \
)


# 1. On lit tous les fichiers
geneBank = myOrthos.GeneBank(noms_fichiers[0], options["speciesUsed"])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], False)


def translateGenome(geneBank, nom, genesAnc):

	# On construit les couleurs
	newGen = {}
	genome = geneBank.dicEspeces[nom]
	for c in genome.lstChr:

		newChr = []
		for gene in genome.lstGenes[c]:
			
			for g in gene.names:
				if g in genesAnc.dicGenes:
					newChr.append( genesAnc.dicGenes[g][1] )
					break
			else:
				newChr.append(nom)
		newGen[c] = newChr
	return newGen


def useDiags(genome1, genome2, nom):

	for c1 in genome1:
		for c2 in genome2:
			if not options["keepGaps"]:
				s1 = set(genome1[c1])
				s2 = set(genome2[c2])
				t1 = [x for x in genome1[c1] if x in s2]
				t2 = [x for x in genome2[c2] if x in s1]
			else:
				t1 = genome1[c1]
				t2 = genome2[c2]
			diags = myMaths.extractDiags(t1, t2, options["fusionThreshold"])
			for d in diags:
				print "%s\t%s" % (nom, d)
	print >> sys.stderr, nom, "OK"


genomes = {}
s = options["speciesUsed"]
for e in s:
	genomes[e] = translateGenome(geneBank, e, genesAnc)

for (i,j) in myTools.myMatrixIterator(len(s), len(s), myTools.myMatrixIterator.StrictUpperMatrix):
	useDiags(genomes[s[i]], genomes[s[j]], s[i]+s[j])

