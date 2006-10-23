#! /usr/bin/python2.4

__doc__ = """
Dessine la matrice des genes orthologues entre deux genomes.
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
import utils.myPsOutput


########
# MAIN #
########

# Arguments
# TODO Couleurs
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["GenomeADessiner", "GenomeReference"], \
	[("taillePoint", float, -1), ("useColor", str, "black"), ("orthologuesList", str, "")], \
	__doc__
)

# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["GenomeADessiner"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["GenomeReference"])
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.AncestralGenome(options["orthologuesList"], False)
else:
	genesAnc = genome2
try:
	colors = utils.myGenomes.loadGenome(options["useColor"])
except Exception:
	colors = options["useColor"]


# Initialisations
nb1 = sum([len(genome1.lstGenes[x]) for x in genome1.lstGenes])
nb2 = sum([len(genome2.lstGenes[x]) for x in genome2.lstGenes])
if options["taillePoint"] < 0:
	dp = 19. / float(nb1)
else:
	dp = options["taillePoint"]


# On ecrit l'entete du PostScipt
utils.myPsOutput.printPsHeader()
utils.myPsOutput.initColor()


def prepareGenome(genome, nb, func):
	i = 0
	y = 1.
	lstNum = {}
	func(y)
	for c in genome.lstChr:
		y += (19. * len(genome.lstGenes[c])) / float(nb)
		func(y)
		for gene in genome.lstGenes[c]:

			for g in gene.names:
				lstNum[g] = i
			i += 1
	return lstNum


# On affiche la grille et on associe "nom de gene" <-> "position sur la grille"
print >> sys.stderr, "Tri des genomes ",
lstNum1 = prepareGenome(genome1, nb1, lambda y: utils.myPsOutput.drawLine(1, y, 19, 0, "black"))
sys.stderr.write(".")
lstNum2 = prepareGenome(genome2, nb2, lambda x: utils.myPsOutput.drawLine(x, 1, 0, 19, "black"))
print >> sys.stderr, ". OK"



print >> sys.stderr, "Affichage des points ",
for c in genome1.lstChr:
	for gene in genome1.lstGenes[c]:
		for g in gene.names:
			x = lstNum1[g]
			if g in genome2.dicGenes:
				(cc,ii) = genome2.dicGenes[g]
				gg = genome2.lstGenes[cc][ii].names
			elif options["orthologuesList"] != "":
				if g in genesAnc.dicGenes:
					(cc,ii) = genesAnc.dicGenes[g]
					gg = genesAnc.lstGenes[cc][ii].names
				else:
					continue
			else:
				continue
			
			for gt in gg:
				if gt in lstNum2:
					y = lstNum2[gt]
					break
			else:
				continue
			
			if type(colors) == str:
				cc = colors
			else:
				for gt in [g]+list(gg):
					if gt in colors.dicGenes:
						cc = colors.dicGenes[gt][0]
						break
				else:
					continue
			xx = 1 + (y*19.)/float(nb2) - dp/2.
			yy = 1 + (x*19.)/float(nb1) - dp/2.
			cc = utils.myPsOutput.getColor(str(cc), "black")
			utils.myPsOutput.drawBox( xx, yy, dp, dp, cc, cc)
			break

utils.myPsOutput.printPsFooter()
print >> sys.stderr, " OK"


