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
	[("KX_HDMC_especes.list",file), ("refGenome",file), ("refSpecies_common",str), ("refSpecies_latin",str), ("diagsFile",file)], \
	[("windowSize",int,10), ("flankLength",int,100)], \
	__doc__ \
)

# Chargement du genome
genome = utils.myGenomes.Genome(arguments["refGenome"])

nbCNE = {}
nbRegions = {}
for i in xrange(-arguments["flankLength"], arguments["flankLength"]+1):
	nbCNE[i] = 0
	nbRegions[i] = 0


posCNE = collections.defaultdict(list)
# Positionnement des CNEs
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	
	if line[0] == arguments["refSpecies_common"]:
		
		c = utils.myGenomes.commonChrName(line[1])
		posCNE[c].append( line[2] )

for x in posCNE.itervalues():
	x.sort()

print >> sys.stderr, sum([len(x) for x in posCNE.itervalues()]), "CNEs"

def addRegion(c, x, x1, x2, coef):

	i1 = bisect.bisect(posCNE[c], x1)
	i2 = min(bisect.bisect(posCNE[c], x2)+1, len(posCNE[c])-1)
	for i in xrange(i1, i2+1):
		xc = posCNE[c][i]
		pos = (coef * (xc - x)) / (arguments["windowSize"]*1000)
		if abs(pos) <= arguments["flankLength"]:
			nbCNE[pos] += 1
	lp = set([(coef * (xc - x)) / (arguments["windowSize"]*1000) for xc in xrange(x1,x2,arguments["windowSize"]*100)])
	for pos in lp:
		if abs(pos) <= arguments["flankLength"]:
			nbRegions[pos] += 1

ft = utils.myFile.myTSV.reader(arguments["diagsFile"])
for t in ft.csvobject:

	# Les genes de la diagonale
	if arguments["refSpecies_latin"] in t[2]:
		do = t[4]
	else:
		assert arguments["refSpecies_latin"] in t[5]
		do = t[7]

	# Les positions de ces genes
	pos = [genome.dicGenes[g] for g in do.split()]
	l = len(pos)
	assert len(set([c for (c,i) in pos])) == 1
	
	(c,imin) = min(pos)
	(_,imax) = max(pos)
	x1 = genome.lstGenes[c][imin].beginning
	x2 = genome.lstGenes[c][imax].end

	if x1 >= x2:
		continue

	addRegion(c, x1, x1-arguments["flankLength"]*arguments["windowSize"]*1000, (x1+x2)/2, 1)
	
	addRegion(c, x2, (x1+x2)/2, x2+arguments["flankLength"]*arguments["windowSize"]*1000, -1)

ft.file.close()

for x in sorted(nbCNE):
	print x, nbCNE[x]/float(nbRegions[x])


