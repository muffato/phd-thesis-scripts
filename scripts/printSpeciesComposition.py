#! /users/ldog/muffato/python -OO

__doc__ = """
Lit des familles de genes et affiche la composition de chaque famille en terme de nombre de genes par espece.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf"],\
	[("speciesList",str,"")], \
	__doc__ \
)

esp = options["speciesList"].split(',')
#geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], only=esp)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"])
print >> sys.stderr, geneBank.dicEspeces.keys()

for s in sys.stdin:
	c = s.split()
	score = dict([(e,0) for e in geneBank.dicEspeces])
	for g in c:
		if g in geneBank.dicGenes:
			(e,_,_) = geneBank.dicGenes[g]
			score[e] += 1
	#print score
	tt = score.values()
	if max(tt) <= 2 and tt.count(2) < tt.count(1):
	#if max(score.values()) >= 2:
		print s, #"\t".join([score[e] for e in esp])
