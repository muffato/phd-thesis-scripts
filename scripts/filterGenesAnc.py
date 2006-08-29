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
	["genesList.conf"],\
	[], \
	"Lit des familles de genes sur l'entree standard et les filtre" \
)


geneBank = myOrthos.MyGeneBank(noms_fichiers[0])

for l in sys.stdin:
	c = l.split()

	score = dict( [(e,0) for e in geneBank.dicEspeces] )

	for g in c:
		if g not in geneBank.dicGenes:
			break
		(e,_,_) = geneBank.dicGenes[g]
		score[e] += 1
	else:
		if score['H'] <= 2 and score['C'] <= 2 and score['M'] <= 2 and score['I'] <= 1 and score['E'] <= 1:
			print l,
