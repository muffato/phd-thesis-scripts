#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import sys
import utils.myGenomes
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs(["GENOME_ANCESTRAL"], [], "")

genesAnc = utils.myGenomes.loadGenome(noms_fichiers["GENOME_ANCESTRAL"])
restrict = set(xrange(100000))

#restrict = [22,32,8,9]
#restrict = [7,29,17,1]
#restrict = [1,21,6,20]
#restrict = [4,19,13,10]
#restrict = [3,23,28,11]
#restrict = [25,5,23,3]

lst = set()
nb = 0
for s in sys.stdin:
	t = set()
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])
	if len(t) < 3:
		continue
	
	nb += 1
	for (c,i) in t:
		if type(c) != int:
			continue
		for (k,j) in t:
			if type(k) != int:
				continue
			if c<k or (c==k and i<j):
				if c not in restrict and k not in restrict:
					continue
				if (c,i,k,j) not in lst:
					print c,i,k,j
					lst.add( (c,i,k,j) )

print >> sys.stderr, nb, "familles pour", len(lst), "couples"
