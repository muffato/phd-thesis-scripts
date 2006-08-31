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
	
for c in genome.lstChr:
	
	y = 0
	
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
		else:
			for g in x:
				if g in genesAnc.dicGenes:
					col = genesAnc.dicGenes[g][0]
					break
			else:
				if not options["includeGaps"]:
					continue
				col = options["defaultColor"]
		col = str(col)
		if col == last:
			nb += 1
		else:
			if nb != 0:
				print "chr%s" % c, y, y+nb-1, "part%s.%d" % (last, y), last
				y += nb
			last = col
			nb = 1
	if nb != 0:
		
		print "chr%s" % c, y, y+nb-1, "part%s.%d" % (last, y), last

