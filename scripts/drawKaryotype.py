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


def loadGenome(nom):
	
	f = open(nom, 'r')
	c = f.readline().split()
	f.close()
	try:
		x = int(c[1]) + int(c[2]) + int(c[3])
		return myOrthos.EnsemblSpeciesGenes(nom)
	except Exception:
		return myOrthos.AncestralGenome(nom, True)


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference"], \
	[("includeGaps", bool, False), ("defaultColor", str, "black")], \
	"Dessine un genome en coloriant ses genes a partir d'un autre genome reference" \
)

genome = loadGenome(noms_fichiers[0])
genesAnc = loadGenome(noms_fichiers[1])
	
dx = (19.*3.)/(5.*len(genome.lstChr)-2.)
dy = 60.

# On ecrit le PostScipt
myPsOutput.printPsHeader()
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
				if not options["includeGaps"]:
					continue
				col = options["defaultColor"]
			else:
				col = genesAnc.dicGenes[x[3]][0]
				col = myPsOutput.getColor(str(col), options["defaultColor"])
		else:
			for g in x:
				if g in genesAnc.dicGenes:
					col = genesAnc.dicGenes[g][0]
					col = myPsOutput.getColor(str(col), options["defaultColor"])
					break
			else:
				if not options["includeGaps"]:
					continue
				col = options["defaultColor"]
		if col == last:
			nb += 1
		else:
			if nb != 0:
				myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
				y += nb/dy
			last = str(col)
			nb = 1
	if nb != 0:
		myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
	xx += (5.*dx)/3.

myPsOutput.printPsFooter()
