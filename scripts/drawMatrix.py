#!/usr/bin/python2.4

# Entree : fichier d'orthologues a 12 champs

# INITIALISATION

# Librairies
import string
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

# Arguments
(noms_fichiers, options) = myTools.checkArgs(["genes_list.conf", "GENES_ANC"], [("taillePoint", float, -1), ("useColors", str, ""), ("espece1", str, "H"), ("espece2", str, "P"), ("ancestralGenome", bool, False)], "")

geneBank = myOrthos.MyGeneBank(noms_fichiers[0], [options["espece1"], options["espece2"]])
if len(geneBank.dicEspeces) != 2:
	print >> sys.stderr, "Can't retrieve %s and %s in %s" % (options["espece1"], options["espece2"], noms_fichiers[0])
	sys.exit(1)
genome1 = geneBank.dicEspeces[options["espece1"]]
genome2 = geneBank.dicEspeces[options["espece2"]]

genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], options["ancestralGenome"])
#dicColors = myOrthos.readDicColors(options["useColors"])

# L'entete pour le dessin de la grille
print "# DEF 1 1 19 19 %f %d %d 0.00001" % (options["taillePoint"], len(genome1.dicGenes), len(genome2.dicGenes))

# Preparation du genome1
i = 0
lstNum1 = {}
s = "# LH 0"
for c in genome1.lstChr:
	s += " " + str(len(genome1.lstGenes[c]))
	for (_,_,_,g) in genome1.lstGenes[c]:
		lstNum1[g] = i
		i += 1
print s

# Preparation du genome2
i = 0
lstNum2 = {}
s = "# LV 0"
for c in genome2.lstChr:
	s += " " + str(len(genome2.lstGenes[c]))
	for (_,_,_,g) in genome2.lstGenes[c]:
		lstNum2[g] = i
		i += 1
print s

# On affiche les genes
for c in genesAnc.lstGenes:
	for g in genesAnc.lstGenes[c]:
		g1 = [x for x in g if x in genome1.dicGenes]
		g2 = set([x for x in g if x in genome2.dicGenes])
		for x in g1:
			for y in g2.difference([x]):
				if options["useColors"] == "":
					print "POINT", lstNum2[y], lstNum1[x]
				else:
					print "POINT", lstNum2[y], lstNum1[x], c
	
sys.exit(0)

# Impression des resultats
for i in range(ortho.nbGenes):

	if options["useColors"] == "":
		print "POINT", ortho.indGenes2[i], ortho.indGenes1[i]
	else:
		# On recupere la couleur
		c = ""
		if options["useColors"] == "1":
			c = ortho.tabGenes1[ortho.indGenes1[i]][0]
		elif options["useColors"] == "2":
			c = ortho.tabGenes2[ortho.indGenes2[i]][0]
		elif options["useColors"] == "3" or options["useColors"] == "4":
			KX = ortho.tabGenes1[ortho.indGenes1[i]][0]
			KY = ortho.tabGenes2[ortho.indGenes2[i]][0]
			c = color[KX][KY]
		elif options["useColors"] != "":
			if ortho.tabGenes1[ortho.indGenes1[i]][3] in dicColors:
				c = dicColors[ortho.tabGenes1[ortho.indGenes1[i]][3]][0]
			elif ortho.tabGenes2[ortho.indGenes2[i]][3] in dicColors:
				c = dicColors[ortho.tabGenes2[ortho.indGenes2[i]][3]][0]
	
		if c != "":
			print "POINT", ortho.indGenes2[i], ortho.indGenes1[i], c


