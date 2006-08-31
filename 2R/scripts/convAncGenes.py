#! /usr/bin/python2.4


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths

(noms_fichiers, options) = myTools.checkArgs(["GENOME_ANCESTRAL"], [], "")

genesAnc = myOrthos.AncestralGenome(noms_fichiers[0], True)

lst = set([])
nb = 0
for s in sys.stdin:
	t = set([])
	#print "*"
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])
	if len(t) >= 2:
		nb += 1

	if 'x' not  in [c for (c,i) in t]:
		continue
	
	for (c,i) in t:
		#print c,
		#continue
		for (k,j) in t:
			if c<k or (c==k and i<j):
				if (c,i,k,j) not in lst:
					print c,i,k,j
					lst.add( (c,i,k,j) )
	#print

print >> sys.stderr, nb, "familles pour", len(lst), "couples conservees"
