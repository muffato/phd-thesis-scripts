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
	[("windowSize",int,10)], \
	__doc__ \
)

# Chargement du genome
genome = utils.myGenomes.Genome(arguments["refGenome"])

tot = collections.defaultdict(float)

# Chargement des diagonales
diags = collections.defaultdict(list)
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
	#assert abs(pos[-1][1]-pos[0][1]) == (l-1)
	(c,imin) = min(pos)
	(_,imax) = max(pos)
	#assert max(pos) == (c,imin+l-1)

	#x1 = genome.lstGenes[c][imin].beginning
	#x2 = genome.lstGenes[c][imin+l-1].end
	#lim = 00000
	#if x2 > x1+2*000000+1:
	#	diags[c].append( (x1+lim,x2-lim) )
	
	#if l >= 15:
	diags[c].append( (genome.lstGenes[c][imin].beginning, genome.lstGenes[c][imax].end) )

for x in diags.itervalues():
	x.sort()

def iterSlices(l):
	while l >= 1000*arguments["windowSize"]:
		yield 1
		l -= 1000*arguments["windowSize"]
	yield float(l)/(1000.*arguments["windowSize"])


# Pour la normalisation
for c in diags:
	# L'interieur des diagonales
	for (x1,x2) in diags[c]:
		if x1 <= x2:
			l = (x2-x1+1)/2
			for (i,x) in enumerate(iterSlices(l)):
				tot[i] += 2*x
		else:
			print >> sys.stderr, "!"
	# Entre des diagonales
	for ((_,x1),(x2,_)) in utils.myTools.myIterator.slidingTuple( diags[c] ):
		if x1 <= x2:
			l = (x2-x1+1)/2
			for (i,x) in enumerate(iterSlices(l)):
				tot[-1-i] += 2*x
		else:
			print >> sys.stderr, "!!"
	# Les deux extremites des chromosomes
	for (i,x) in enumerate(iterSlices(diags[c][0][0])):
		tot[-1-i] += x
	for (i,x) in enumerate(iterSlices(10*1000*arguments["windowSize"])):
		tot[-1-i] += x


res = collections.defaultdict(int)
# Positionnement des CNEs
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	
	if line[0] == arguments["refSpecies_common"]:
		
		c = utils.myGenomes.commonChrName(line[1])
		i = bisect.bisect_left(diags[c], (line[2],line[3]))
		
		if (i >= 1) and (line[2] < diags[c][i-1][1]):
			# dans diagonale precedente
			choix1 = diags[c][i-1][0]
			choix2 = diags[c][i-1][1]
			into = True
		else:
			# entre diagonales
			choix1 = diags[c][i-1][1] if i >= 1 else 0
			choix2 = diags[c][i][0] if i < len(diags[c]) else sys.maxint
			into = False

		assert choix1 < line[2] < line[3], (choix1,line[2], line[3],choix2)

		nb1 = (line[2]-choix1) / (1000*arguments["windowSize"])
		nb2 = (choix2-line[3]) / (1000*arguments["windowSize"])
		if nb2 == -1:
			nb2 = 0
		assert min(nb1,nb2) >= 0
		
		add = 1

		if into:
			res[min(nb1,nb2)] += add
		else:
			res[-1-min(nb1,nb2)] += add

print >> sys.stderr, sum(tot.values())
for x in sorted(res):
	print utils.myFile.myTSV.printLine( [x, res[x], tot[x], res[x]/tot[x]] )

