#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit la matrice des p-values et extrait les clusters de fortement lies
"""

import sys
import bisect
import utils.myTools

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["pvalueFile"], \
	[("cutoff",float,5)], \
	__doc__
)

f = utils.myTools.myOpenFile(noms_fichiers["pvalueFile"], "r")
pvalues = utils.myTools.defaultdict(dict)
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
			if w >= options["cutoff"]:
				print '%s -- %s [weight="%f"]' % (c1,c2,w)
print "}"

