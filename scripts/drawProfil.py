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
	["GenomeADessiner", "GenomeReference"], \
	[], \
	"Dessine un genome en coloriant ses genes a partir d'un autre genome reference" \
)


if "genes" in noms_fichiers[0]:
	genome = myOrthos.EnsemblSpeciesGenes(noms_fichiers[0])
else:
	genome = myOrthos.AncestralGenome(noms_fichiers[0], True)
	
if "genes" in noms_fichiers[1]:
	genesAnc = myOrthos.EnsemblSpeciesGenes(noms_fichiers[1])
else:
	genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], True)
	
dx = (19.*3.)/(5.*len(genome.lstChr)-2.)
dy = 60.

# On ecrit le PostScipt
myPsOutput.printPsHeader()
myPsOutput.loadColorTable()
myPsOutput.initColor()


xx = 1
y0 = 1.

for c in genome.lstChr:
	
	myPsOutput.drawText(xx, y0, str(c), "black")
	y = y0 + 1
	
	last = ""
	nb = 0
	for x in genome.lstGenes[c]:
		if type(x) == tuple:
			if x[3] not in genesAnc.dicGenes:
				continue
			col = genesAnc.dicGenes[x[3]][0]
		else:
			for g in x:
				if g in genesAnc.dicGenes:
					(col,_) = genesAnc.dicGenes[g]
					break
			else:
				continue
		
		if col == last:
			nb += 1
		else:
			if nb != 0:
				c = myPsOutput.getColor(last, 0)
				myPsOutput.drawBox(xx, y, dx, nb/dy, c, c)
				y += nb/dy
			last = str(col)
			nb = 1
	if nb != 0:
		c = myPsOutput.getColor(last, 0)
		myPsOutput.drawBox(xx, y, dx, nb/dy, c, c)
	xx += (5.*dx)/3.

myPsOutput.printPsFooter()
