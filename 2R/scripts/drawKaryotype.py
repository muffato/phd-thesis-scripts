#! /usr/bin/python2.4

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
import myPsOutput

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference", "orthologues"], \
	[("includeGaps", bool, False), ("defaultColor", str, "black")], \
	"Dessine un genome en coloriant ses genes a partir d'un autre genome reference" \
)

genome = myOrthos.EnsemblGenome(noms_fichiers[0])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], True)
orthologues = myOrthos.AncestralGenome(noms_fichiers[2], False)

# On ecrit le PostScipt
myPsOutput.printPsHeader()
myPsOutput.initColor()

# On calcule les couleurs
res = {}
for c in genome.lstChr:
	
	res[c] = []
	
	for x in genome.lstGenes[c]:
		if x[3] not in orthologues.dicGenes:
			if not options["includeGaps"]:
				continue
			col = options["defaultColor"]
		else:
			(cc,ii) = orthologues.dicGenes[x[3]]
			t = orthologues.lstGenes[cc][ii]
			for s in t:
				if s in genesAnc.dicGenes:
					col = genesAnc.dicGenes[s][0]
					col = myPsOutput.getColor(str(col), options["defaultColor"])
					break
			else:
				if not options["includeGaps"]:
					continue
				col = options["defaultColor"]
		res[c].append(col)

# On dessine 
dx = (19.*3.)/(5.*len(genome.lstChr)-2.)
dy = float(max([len(x) for x in res.values()])) / 26.
xx = 1
y0 = 1.

for c in genome.lstChr:
	
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
