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

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths
import myPsOutput

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference"], \
	[("includeGaps", bool, False), ("defaultColor", str, "black"), ("orthologuesList", str, "")], \
	__doc__ \
)

genome1 = myOrthos.loadGenome(noms_fichiers[0])
genome2 = myOrthos.loadGenome(noms_fichiers[1])
if options["orthologuesList"] != "":
	genesAnc = myOrthos.AncestralGenome(options["orthologuesList"], False)

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
			tg = myMaths.flatten([genesAnc.lstGenes[cc][ii].names for (cc,ii) in [genesAnc.dicGenes[g] for g in tg if g in genesAnc.dicGenes]])
		
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
sys.exit(0)
# On dessine
dx = (19.*3.)/(5.*len(genome1.lstChr)-2.)
dy = float(max([len(x) for x in res.values()])) / 26.
xx = 1
y0 = 1.


for c in genome1.lstChr:
	
	myPsOutput.drawText(xx, y0, str(c), "black")
	y = y0 + 1
	
	last = ""
	nb = 0
	for col in res[c]:
		if col == last:
			nb += 1
		else:
			if nb != 0:
				myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
				y += nb/dy
			last = col
			nb = 1
	if nb != 0:
		myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
	xx += (5.*dx)/3.

myPsOutput.printPsFooter()

