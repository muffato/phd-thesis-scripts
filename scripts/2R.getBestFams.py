#!/usr/bin/env python2

__doc__ = """
	Lit la matrice des p-values et extrait les clusters de fortement lies
"""

import sys
import bisect
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs( \
	[("pvalueFile",file)], \
	[("clustersSizes",str,"2,3,4,5,6,7"), ("cutoff",float,5)], \
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

res = []
for n in arguments["clustersSizes"].split(","):
	n = int(n)
	for t in utils.myTools.myIterator.buildSubsets(chrom, n):
		t = tuple(t)
		s = 0.
		nb = 0
		for (x1,x2) in utils.myTools.myIterator.tupleOnStrictUpperList(t):
			s += pvalues[x1][x2]
			nb += 1
		s /= nb
		if s >= arguments["cutoff"]:
			bisect.insort(res, (s,nb,t))
	print >> sys.stderr, n, "ok"
print >> sys.stderr, len(res)

for x in res:
	print utils.myFile.myTSV.printLine(x)

