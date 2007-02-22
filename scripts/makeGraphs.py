#! /users/ldog/muffato/python -OO

__doc__ = """
Lit sur l'entree standard des familles d'objets et les regroupe.
Il ne s'agit que d'une interface pour l'objet myCombinator.
"""

import sys
import math
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs(["familiesFile"], [], __doc__)

# MAIN #

# On scanne toutes les lignes
f = utils.myTools.myOpenFile(noms_fichiers["familiesFile"], 'r')
dic = {}
n = 0
for l in f:
	t = [int(x) for x in l.split()]
	for i in xrange(len(t)):
		dic[t[i]] = (n,i)
	n += 1
f.close()

n = 0
lastG = -1
f = None
for l in sys.stdin:
	
	t = l.split()

	a = int(t[1])
	b = int(t[2])
	c = float(t[11])
	if c == 0:
		c = 1000000
	else:
		c = -math.log10(c)

	(g,ia) = dic[a]
	(_,ib) = dic[b]

	if g != lastG:
		if f != None:
			f.close()
		f = open("graph.%d" % g, 'a')
		lastG = g

	print >> f, ia, ib, c

	n += 1
	if n % 100000 == 0:
		print >> sys.stderr, n
