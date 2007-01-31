#! /users/ldog/muffato/python -OO

__doc__ = """
Lit des familles de genes sur l'entree standard.
Calcule le nombre de genes de chaque espece et agit en consequence
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import os
import utils.myGenomes
import utils.myTools


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf"],\
	[("breakWhenFamilyNotComplete",bool,False), ("speciesList",str,"")], \
	__doc__ \
)

esp = options["speciesList"].split(',')
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], only=esp)


for l in sys.stdin:
	c = l.split()

	score = dict( [(e,0) for e in geneBank.dicEspeces] )

	for g in c:
		if g not in geneBank.dicGenes:
			if options["breakWhenFamilyNotComplete"]:
				break
			else:
				continue
		(e,_,_) = geneBank.dicGenes[g]
		score[e] += 1
	else:
		#t = score.values()
		t = [score[x] for x in geneBank.lstEspecesNonDup]
		tt = [score[x] for x in geneBank.lstEspecesDup]
		#tt = [score[x] for x in ]
		#if score['H'] == 1 and score['C'] == 1:
		#if max(t) > 0 and max(tt) > 0:
		#if max(t) == 2 and max(tt) == 1:
		#if sum(t) < 13 and sum(tt) < 6:
		#if len(t) == 1 and 1 in t:
		if max(t) <= 1: # and max(tt) <= 2:
		#if max(t) == 1 and min(t) == 1:
			print l,

