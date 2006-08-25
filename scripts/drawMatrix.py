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
	["genesList.conf", "genesAncestraux.list"],\
	[("taillePoint", float, -1), ("useColors", bool, True), ("espece1", str, ""), ("espece2", str, ""), ("ancestralGenome", bool, False)], \
	"Dessine une matrice d'orthologie a partir de la liste des genes orthologues et des genomes des especes" \
)


geneBank = myOrthos.MyGeneBank(noms_fichiers[0], [options["espece1"], options["espece2"]])

if len(geneBank.dicEspeces) == 0 or (len(geneBank.dicEspeces) == 1 and not options["ancestralGenome"]):
	print >> sys.stderr, "Can't retrieve -%(espece1)s- and -%(espece2)s- in the species list" % options
	sys.exit(1)

genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], options["ancestralGenome"])
if options["espece1"] != "":
	genome1 = geneBank.dicEspeces[options["espece1"]]
else:
	genome1 = genesAnc
nb1 = sum([len(genome1.lstGenes[x]) for x in genome1.lstGenes])

if options["espece2"] != "":
	genome2 = geneBank.dicEspeces[options["espece2"]]
else:
	genome2 = genesAnc
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
		if genome == genesAnc:
			for s in genome.lstGenes[c]:
				for g in s:
					lstNum[g] = i
				i += 1
		else:
			for (_,_,_,g) in genome.lstGenes[c]:
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
for c in genesAnc.lstGenes:
	for g in genesAnc.lstGenes[c]:
		if genome1 == genesAnc:
			g1 = [g[0]]
		else:
			g1 = [x for x in g if x in genome1.dicGenes]
		if genome2 == genesAnc:
			g2 = [g[0]]
		else:
			g2 = [x for x in g if x in genome2.dicGenes]
		for x in g1:
			for y in g2:
				xx = 1 + (lstNum2[y]*19.)/float(nb2) - dp/2.
				yy = 1 + (lstNum1[x]*19.)/float(nb1) - dp/2.

				if options["useColors"]:
					cc = myPsOutput.getColor(c, "black")
				else:
					cc = "black"
				myPsOutput.drawBox( xx, yy, dp, dp, cc, cc)
	sys.stderr.write(".")

myPsOutput.printPsFooter()
print >> sys.stderr, " OK"

