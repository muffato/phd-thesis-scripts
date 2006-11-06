#! /usr/bin/python2.4

__doc__ = """
Dessine un genome en coloriant ses genes a partir d'un autre genome reference.
"""

##################
# INITIALISATION #
##################

# Librairies
import string
import sys
import os
import utils.myOrthos
import utils.myTools
import utils.myMaths
import utils.myPsOutput

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference"], \
	[("includeGaps", bool, False), ("defaultColor", str, "black"), ("orthologuesList", str, "")], \
	__doc__ \
)

genome1 = utils.myOrthos.loadGenome(noms_fichiers[0])
genome2 = utils.myOrthos.loadGenome(noms_fichiers[1])
if options["orthologuesList"] != "":
	genesAnc = utils.myOrthos.AncestralGenome(options["orthologuesList"], False)

# On ecrit le PostScipt
#myPsOutput.printPsHeader()
#myPsOutput.initColor()

# On construit les couleurs
res = {}
somme = 0
total = 0
for c in genome1.lstChr:

	res[c] = []
	for gene in genome1.lstGenes[c]:
		
		tg = gene.names
		if options["orthologuesList"] != "":
			tg = utils.myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in tg if g in genesAnc.dicGenes]])
		
		for g in tg:
			if g in genome2.dicGenes:
				col = genome2.dicGenes[g][0]
				#col = myPsOutput.getColor(str(col), options["defaultColor"])
				break
		else:
			#if not options["includeGaps"]:
			#	continue
			col = options["defaultColor"]
		res[c].append(col)
		total += 1
		if col == c:
			somme += 1

print >> sys.stderr, somme, total

