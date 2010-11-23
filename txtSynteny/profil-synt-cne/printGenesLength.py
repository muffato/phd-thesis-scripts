#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Lit les CNE et transforme les positions genomiques en positions indexees
"""

import sys
import bisect
import operator
import itertools
import collections

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( \
	[("intervalFile",file)], \
	[], \
	__doc__ \
)



# Chargement du genome
print >> sys.stderr, "Loading intervals ...",
genome = collections.defaultdict(list)
f = utils.myFile.openFile(arguments["intervalFile"], "r")
for l in f:
	l = eval(l)
	genome[l[0]].append(l[1:])
f.close()

dicGene2Pos = collections.defaultdict(list)
for c in genome:
	#assert genome[c] == sorted(genome[c])
	#for (x,y) in utils.myTools.myIterator.slidingTuple(genome[c]):
	#	assert y[0] == x[1]+1
	for (i,x) in enumerate(genome[c]):
		for g in x[3]:
			dicGene2Pos[g].append( (c,i) )
print >> sys.stderr, "OK"


# Chargement des diagonales
print >> sys.stderr, "Loading pseudo-synteny blocks ...",
for x in dicGene2Pos:

	if isinstance(x, frozenset):
		continue

	(cmin,imin) = min(dicGene2Pos[x])
	(cmax,imax) = max(dicGene2Pos[x])
	assert cmin == cmax
	c = cmin
	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	print xmax-xmin+1

print >> sys.stderr, "OK"

