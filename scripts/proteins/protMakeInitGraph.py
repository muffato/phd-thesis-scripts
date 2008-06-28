#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit les scores d'alignement d'Ensembl et affiche cree le graphe correspondant (+uniq) et imprime les composantes connexes
"""

import sys
import math
import utils.myTools


arguments = utils.myTools.checkArgs( [], [], __doc__)


comb = utils.myTools.myCombinator([])

for l in sys.stdin:
	
	t = l.split()

	a = int(t[1])
	b = int(t[2])
	
	if a == b:
		continue
	
	c = float(t[11])
	if c == 0:
		c = 1000000
	else:
		c = -math.log10(c)

	print a, b, c
	
	comb.addLink([a,b])

for g in comb:
	print >> sys.stderr, utils.myTools.printLine(g, " ")

