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
import utils.myGenomes
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

genome1 = utils.myGenomes.loadGenome(noms_fichiers["GenomeADessiner"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["GenomeReference"])
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.AncestralGenome(options["orthologuesList"], False, False)

# On ecrit le PostScipt
utils.myPsOutput.printPsHeader()
utils.myPsOutput.initColor()

# On construit les couleurs
res = {}
for c in genome1.lstChr:

	res[c] = []
	for gene in genome1.lstGenes[c]:
		
		tg = gene.names
		if options["orthologuesList"] != "":
			tg = utils.myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in tg if g in genesAnc.dicGenes]])
		
		for g in tg:
			if g in genome2.dicGenes:
				col = genome2.dicGenes[g][0]
				col = utils.myPsOutput.getColor(str(col), options["defaultColor"])
				break
		else:
			if not options["includeGaps"]:
				continue
			col = options["defaultColor"]
			#print >> sys.stderr, gene.names
		res[c].append(col)

# On dessine
dx = (19.*3.)/(5.*len(genome1.lstChr)-2.)
dy = float(max([len(x) for x in res.values()])) / 26.
xx = 1
y0 = 1.

for c in genome1.lstChr:
	
	utils.myPsOutput.drawText(xx, y0, str(c), "black")
	y = y0 + 1
	
	last = ""
	nb = 0
	for col in res[c]:
		if col == last:
			nb += 1
		else:
			if nb != 0:
				utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
				y += nb/dy
			last = col
			nb = 1
	if nb != 0:
		utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
	xx += (5.*dx)/3.

utils.myPsOutput.printPsFooter()

