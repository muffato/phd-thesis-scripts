#! /usr/bin/python2.4

##################
# INITIALISATION #
##################

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["ancestralGenome"],\
	[], \
	"Lit des familles de genes sur l'entree standard et les filtre" \
)


geneBank = myOrthos.AncestralGenome(noms_fichiers[0], True)

for l in sys.stdin:
	c = l.split()

	if c[0] in geneBank.dicGenes:
		print geneBank.dicGenes[c[0]][0]
