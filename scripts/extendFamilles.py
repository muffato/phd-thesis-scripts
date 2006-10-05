#! /usr/bin/python2.4

__doc__ = """
Lit des familles de genes sur l'entree standard.
Rajoute a chaque famille les genes associes dans le genome ancestral passe en argument.
"""


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools

(noms_fichiers, options) = myTools.checkArgs(["genesAncestraux"], [], __doc__)

genesAnc = myOrthos.loadGenome(noms_fichiers[0])

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
