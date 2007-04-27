#! /users/ldog/muffato/python -OO

##################
# INITIALISATION #
##################

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myGenomes
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


geneBank = myGenomes.loadGenome(noms_fichiers["ancestralGenome"])

for l in sys.stdin:
	c = l.split()

	if c[0] in geneBank.dicGenes:
		print geneBank.dicGenes[c[0]][0]
