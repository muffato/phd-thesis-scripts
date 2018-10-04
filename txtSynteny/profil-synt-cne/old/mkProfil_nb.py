#!/usr/bin/env python2

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
	[("flankLength",int,10)], \
	__doc__ \
)

# Chargement du genome
genome = utils.myGenomes.Genome(arguments["refGenome"])

# Initialisation des comptages
nbCNE = {}
longCNE = {}
longGenome = {}
for (c,x) in genome.lstGenes.iteritems():
	d = {}
	dl = {}
	lx = len(x)
	for i in xrange(lx):
		d[i] = 0
		dl[i] = x[i].end - x[i].beginning
		d[(i-1,i)] = 0
		if i > 0:
			dl[(i-1,i)] = x[i].beginning - x[i-1].end
			#if dl[(i-1,i)] < 0:
			#	print dl[(i-1,i)], x[i-1], x[i]
		else:
			dl[(-1,0)] = x[0].beginning
	d[(lx-1,lx)] = 0
	dl[(lx-1,lx)] = 0 # devrait etre la distance jusqu'a la fin du chromosome

	for j in xrange(arguments["flankLength"]+1):
		d[-1-j] = 0
		d[(-2-j,-1-j)] = 0
		d[lx+j] = 0
		d[(lx+j,lx+1+j)] = 0
		dl[-1-j] = 0
		dl[(-2-j,-1-j)] = 0
		dl[lx+j] = 0
		dl[(lx+j,lx+1+j)] = 0
	nbCNE[c] = d.copy()
	longCNE[c] = d.copy()
	longGenome[c] = dl



# Positionnement des CNEs
for line in utils.myFile.myTSV.readTabular(arguments["KX_HDMC_especes.list"], [str,str,int,int,str,str,int,str]):
	
	if line[0] == arguments["refSpecies_common"]:
		
		c = utils.myGenomes.commonChrName(line[1])

		i = bisect.bisect_left(genome.lstGenes[c], (line[2],line[3]))
		i = bisect.bisect_left(genome.lstGenes[c], utils.myGenomes.Gene(c, line[2], line[3], int(line[4]+"1"), line[5]))

		if (i >= 1) and (genome.lstGenes[c][i-1].end > line[2]):
			# intronique
			pos = i-1
		else:
			# inter genique
			pos = (i-1,i)
		
		nbCNE[c][pos] += 1
		longCNE[c][pos] += line[3]-line[2]+1


res = collections.defaultdict(list)
def printInterval(j, c, i):
	print "!", utils.myFile.myTSV.printLine( [(j,j+1), nbCNE[c][(i,i+1)], longCNE[c][(i,i+1)], longGenome[c][(i,i+1)], genome.lstGenes[c][i] if i < len(genome.lstGenes[c]) and i >= 0 else None, genome.lstGenes[c][i+1] if i+1 < len(genome.lstGenes[c]) and i >= 0 else None] )
	res[(j,j+1)].append( (nbCNE[c][(i,i+1)], longCNE[c][(i,i+1)], longGenome[c][(i,i+1)], (genome.lstGenes[c][i].strand*genome.lstGenes[c][i+1].strand) if i+1 < len(genome.lstGenes[c]) and i >= 0 else None  ) )

def printGene(j, c, i, mul):
	print "!", utils.myFile.myTSV.printLine( [j, nbCNE[c][i], longCNE[c][i], longGenome[c][i], genome.lstGenes[c][i] if i < len(genome.lstGenes[c]) and i >= 0 else None] )
	res[j].append( (nbCNE[c][i], longCNE[c][i], longGenome[c][i], genome.lstGenes[c][i].strand*mul if i < len(genome.lstGenes[c]) and i >= 0 else None) )

todoI = collections.defaultdict(list)
seenI = set()
todoG = collections.defaultdict(list)
seenG = set()
posBlocks = collections.defaultdict(set)
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
	l = max(pos)[1]-imin+1
	#assert max(pos) == (c,imin+l-1)

	# Parcours des intervalles
	for j in xrange(-arguments["flankLength"],l+1+arguments["flankLength"]):
		i = imin + j - 1
		j = min(j, l-j)
		if j < 0:
			todoI[(c,i)].append(j)
		elif (c,i) not in seenI:
			seenI.add( (c,i) )
			printInterval(j, c, i)
		#else:
		#	assert j in [0,1]

	# Parcours des genes
	for j in xrange(-arguments["flankLength"],l+2+arguments["flankLength"]):
		i = imin + j - 1
		lj = l+1-j
		j = min(j, l+1-j)
		mul = -1 if j == lj else 1
		if j <= 0:
			todoG[(c,i)].append((j,mul))
		elif (c,i) not in seenG:
			seenG.add( (c,i) )
			printGene(j, c, i, mul)
		#else:
		#	assert j == 1

ft.file.close()


# Intervales flanquants non synteniques
for ((c,i),l) in todoI.iteritems():
	if (c,i) not in seenI:
		printInterval(max(l), c, i)

# Genes flanquants non synteniques
for ((c,i),l) in todoG.iteritems():
	if (c,i) not in seenG:
		(l,mul) = max(l)
		printGene(l, c, i, mul)

toprint = range(-500, 500) + [(i,i+1) for i in xrange(-500, 500)]

for key in toprint:

	if key not in res:
		continue
	lst = res[key]
	
	#lst = sorted(res[key])[:-len(lst)/10]
	
	lst = [x for x in lst if x[2]>0 ]
	if len(lst) == 0:
		continue

	totNbCNE = sum(x[0] for x in lst)
	totLongCNE = sum(x[1] for x in lst)
	totLongGenome = sum(x[2] for x in lst)
	strands = [x[3] for x in lst].count(1)
	
	if totLongGenome == 0:
		continue

	nb = float(len(lst))
	print utils.myFile.myTSV.printLine([key, totNbCNE/nb, totLongCNE/float(totNbCNE) if totNbCNE != 0 else None, totLongGenome/nb, (1000000.*totNbCNE)/totLongGenome, (100.*totLongCNE)/totLongGenome, 100.*strands/nb, len(lst)])



