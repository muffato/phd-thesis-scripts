#!/usr/bin/env python2

__doc__ = """
	Lit la matrice des p-values et cree le graphe (en fonction du seuil)
"""

import sys
import bisect
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs( \
	[("pvalueFile",file)], \
	[("cutoff",float,5)], \
	__doc__
)

f = utils.myTools.myOpenFile(arguments["pvalueFile"], "r")
pvalues = collections.defaultdict(dict)
for ligne in f:
	t = ligne.split()
	if t[0] == '+':
		chrom = t[1:]
	else:
		for (i,x) in enumerate(t[1:]):
			pvalues[t[0]][chrom[i]] = float(x)
f.close()

print "graph {"

for c1 in pvalues:
	for (c2,w) in pvalues[c1].iteritems():
		if c2 > c1:
			if w >= arguments["cutoff"]:
				print '%s -- %s [weight="%f"]' % (c1,c2,w)
print "}"

