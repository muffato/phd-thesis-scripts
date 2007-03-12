#! /users/ldog/muffato/python -OO

__doc__ = """
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
dejaLus = set([])

espOK = set([3, 4, 13, 15, 22, 25, 26, 28, 31, 36, 37, 38, 39, 42])

for l in sys.stdin:
	
	t = l.split()

	if (int(t[3]) not in espOK) or (int(t[4]) not in espOK):
		continue
		
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

	(g,ia) = dic[a]
	(_,ib) = dic[b]

	if g != lastG:
		if f != None:
			f.close()
		f = open("graph.%d" % g, 'a')
		lastG = g

	print >> f, ia, ib, c

