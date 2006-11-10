#! /users/ldog/muffato/python

__doc__ = """
Lit des familles de genes sur l'entree standard.
Rajoute a chaque famille les genes associes dans le genome ancestral passe en argument.
"""


# INITIALISATION #

# Librairies
import sys
import os
import utils.myGenomes
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs(["genesAncestraux"], [], __doc__)

genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genesAncestraux"])

lst = set([])
nb = 0
for s in sys.stdin:
	t = set([])
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])
	
	for (c,i) in t:
		print " ".join(genesAnc.lstGenes[c][i].names),
	print s,

print >> sys.stderr, nb, "familles pour", len(lst), "couples conservees"
