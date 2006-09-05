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
	
	f = myTools.myOpenFile(nom)
	c = f.readline().split()
	f.close()
	try:
		x = int(c[1]) + int(c[2]) + int(c[3])
		return myOrthos.EnsemblGenome(nom)
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

# On ecrit le PostScipt
myPsOutput.printPsHeader()
myPsOutput.initColor()

# On construit les couleurs
res = {}
for c in genome.lstChr:

	res[c] = []
	
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
		#print c, x, col
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
