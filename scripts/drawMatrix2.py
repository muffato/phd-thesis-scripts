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
	["genesAncestraux1.list", "genesAncestraux2.list"],\
	[("taillePoint", float, -1)], \
	"Dessine une matrice d'orthologie a partir de la liste des genes orthologues et des genomes des especes" \
)


genome1 = myOrthos.AncestralGenome(noms_fichiers[0], True)
genome2 = myOrthos.AncestralGenome(noms_fichiers[1], True)
nb1 = sum([len(genome1.lstGenes[x]) for x in genome1.lstGenes])
nb2 = sum([len(genome2.lstGenes[x]) for x in genome2.lstGenes])


# On ecrit le PostScipt
myPsOutput.printPsHeader()
myPsOutput.initColor()

if options["taillePoint"] < 0:
	dp = 19. / float(nb1)
else:
	dp = options["taillePoint"]


def prepareGenome(genome, nb, func):
	i = 0
	y = 1.
	lstNum = {}
	func(y)
	for c in genome.lstChr:
		y += (19. * len(genome.lstGenes[c])) / float(nb)
		func(y)
		for s in genome.lstGenes[c]:
			for g in s:
				lstNum[g] = i
			i += 1
	return lstNum



# On affiche les genes
print >> sys.stderr, "Tri des genomes ",
lstNum1 = prepareGenome(genome1, nb1, lambda y: myPsOutput.drawLine(1, y, 19, 0, "black"))
sys.stderr.write(".")
lstNum2 = prepareGenome(genome2, nb2, lambda x: myPsOutput.drawLine(x, 1, 0, 19, "black"))
print >> sys.stderr, ". OK"

print >> sys.stderr, "Affichage des points ",
for c in genome1.lstChr:
	for g in genome1.lstGenes[c]:
		g1 = [g[0]]
		g2 = [g[0]]
		#g1 = [x for x in g if x in genome1.dicGenes]
		for x in g1:
			for y in g2:
				xx = 1 + (lstNum2[y]*19.)/float(nb2) - dp/2.
				yy = 1 + (lstNum1[x]*19.)/float(nb1) - dp/2.

				cc = "black"
				myPsOutput.drawBox( xx, yy, dp, dp, cc, cc)
	sys.stderr.write(".")

myPsOutput.printPsFooter()
print >> sys.stderr, " OK"

