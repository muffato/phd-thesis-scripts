#! /users/ldog/muffato/python

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
	[("intervalFile",file), ("refSpecies_latin",str), ("diagsFile",file)], \
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
print >> sys.stderr, "Loading synteny blocks ...",
ft = utils.myFile.myTSV.reader(arguments["diagsFile"])
for t in ft.csvobject:

	# Les genes de la diagonale
	if arguments["refSpecies_latin"] in t[2]:
		c = utils.myGenomes.commonChrName(t[3])
		do = t[4]
	else:
		assert arguments["refSpecies_latin"] in t[5]
		c = utils.myGenomes.commonChrName(t[6])
		do = t[7]

	# Les positions de ces genes
	genes = do.split()
	lg = len(genes)
	ind = dict([(g,i) for (i,g) in enumerate(genes)])
	(cmin,imin) = min([min(dicGene2Pos[g]) for g in genes])
	(cmax,imax) = max([max(dicGene2Pos[g]) for g in genes])
	assert c == cmin
	assert c == cmax

	xmin = genome[c][imin][0]
	xmax = genome[c][imax][1]
	print lg, xmax-xmin+1
ft.file.close()
print >> sys.stderr, "OK"

