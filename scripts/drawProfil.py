#!/usr/bin/python2.4

# Entree : fichier d'orthologues a 12 champs

# INITIALISATION

# Librairies
import string
import sys
import os

sys.path.append(os.environ['HOME'] + "/M2/scripts/utils")
import myOrthos
import myTools

# Arguments
(noms_fichiers, options) = myTools.checkArgs(["data/genes_E", "GENOME_ANCESTR"], [("ancestralColors", bool, True)])


genome = myOrthos.EnsemblSpeciesGenes(noms_fichiers[0])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], True)

if options["ancestralColors"]:
	tabC = genome.lstChr
	tabG = genome.lstGenes
else:
	tabC = genesAnc.lstChr
	tabG = genesAnc.lstGenes

x = (20.*3.)/(5.*(2+len(tabC))-2.)
print ">Z 0.5 3 %f %f 50" % (x, 2.*x/3.)
print


for c in tabC:
	s = str(c)
	last = ""
	nb = 0
	for x in tabG[c]:
		if options["ancestralColors"]:
			if x[3] not in genesAnc.dicGenes:
				continue
			col = genesAnc.dicGenes[x[3]][0]
		else:
			for g in x:
				if g in genome.dicGenes:
					(col,_) = genome.dicGenes[g]
					break
			else:
				continue
		
		if col == last:
			nb += 1
		else:
			if nb != 0:
				s += " " + last + " " + str(nb)
			last = str(col)
			nb = 1
	if nb != 0:
		s += " " + last + " " + str(nb)
	print s


