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


geneBank = myOrthos.GeneBank(noms_fichiers[0])

for l in sys.stdin:
	c = l.split()

	score = dict( [(e,0) for e in geneBank.dicEspeces] )

	for g in c:
		if g not in geneBank.dicGenes:
			continue
			#break
		(e,_,_) = geneBank.dicGenes[g]
		score[e] += 1
	else:
		t = [score[x] for x in 'HMDWOC']
		tt = [score[x] for x in 'TSZF']
		#if score['H'] == 1 and score['C'] == 1:
		#if max(t) > 0 and max(tt) > 0:
		#if max(t) == 2 and max(tt) == 1:
		#if sum(t) < 13 and sum(tt) < 6:
		#if len(t) == 1 and 1 in t:
		if max(t) == 1 and max(tt) == 2:
		#if max(t) == 1 and min(t) == 1:
			print l,
