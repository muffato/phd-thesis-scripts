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
	["genesList.conf"],\
	[], \
	"Lit des familles de genes et affiche la composition de chaque famille en terme de nombre de genes par espece" \
)


geneBank = myOrthos.GeneBank(noms_fichiers[0])

for s in sys.stdin:
	c = s.split()
	score = dict([(e,0) for e in geneBank.dicEspeces])
	for g in c:
		if g in geneBank.dicGenes:
			(e,_,_) = geneBank.dicGenes[g]
			score[e] += 1
	if score['H'] > 1 or score['C'] > 1 or score['M'] > 1 or score['O'] > 1 or score['D'] > 1 or score['W'] > 1:
		print score
	#print score['H']
