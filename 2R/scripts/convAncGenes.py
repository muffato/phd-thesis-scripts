#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myGenomes
import myTools
import myMaths

(noms_fichiers, options) = myTools.checkArgs(["GENOME_ANCESTRAL"], [], "")

genesAnc = myGenomes.loadGenome(noms_fichiers["GENOME_ANCESTRAL"])

lst = set([])
nb = 0
for s in sys.stdin:
	t = set([])
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])

	if len(t) > 0:
		nb += 1

	for (c,i) in t:
		if type(c) != int:
			continue
		for (k,j) in t:
			if type(k) != int:
				continue
			if c<k or (c==k and i<j):
				if (c,i,k,j) not in lst:
					print c,i,k,j
					lst.add( (c,i,k,j) )

print >> sys.stderr, nb, "familles pour", len(lst), "couples"
