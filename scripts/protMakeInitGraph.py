#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit un graphe (aretes) et supprime la redondance et imprime les composantes connexes
"""

import sys
import math
import utils.myTools


(noms_fichiers, options) = utils.myTools.checkArgs( [], [], __doc__)


dejaLus = set([])
comb = utils.myTools.myCombinator([])

for l in sys.stdin:
	
	t = l.split()

	a = int(t[1])
	b = int(t[2])
	
	if a == b:
		continue
	
	elif a > b:

		if (b,a) in dejaLus:
			continue

		g = a
		a = b
		b = g

	else:
		dejaLus.add( (a,b) )
	
	c = float(t[11])
	if c == 0:
		c = 1000000
	else:
		c = -math.log10(c)

	print a, b, c
	
	comb.addLink([a,b])

dejaLus = None
for g in comb:
	print >> sys.stderr, " ".join([str(x) for x in g])

