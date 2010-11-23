#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Cherche des points de rearrangements par syntenie homme/ancetre
"""

import sys
import itertools
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import myDiags

# Arguments
arguments = utils.myTools.checkArgs( \
	[("ancGenes",file), ("ancGenome",file), ("modernGenome",file)], \
	[("minimalLength",int,2)], \
	__doc__ \
)

ancGenes = utils.myGenomes.Genome(arguments["ancGenes"])
ancGenome = utils.myGenomes.Genome(arguments["ancGenome"], ancGenes=ancGenes)
modernGenome = utils.myGenomes.Genome(arguments["modernGenome"])

# Les genes sur des chromosomes / scaffolds / blocs trop petits ne peuvent creer de syntenie utilisable
n = 0
while n != len(ancGenes.lstGenes[None]):
	n = len(ancGenes.lstGenes[None])
	ref = set(ancGenes.dicGenes)
	filterOut = set()
	def filt(genome):
		for (chrom,l) in genome.lstGenes.iteritems():
			c = 0
			for gene in l:
				if len(ref.intersection(gene.names)) > 0:
					c += 1
				if c == arguments["minimalLength"]:
					break
			else:
				for gene in l:
					filterOut.update(gene.names)
	filt(ancGenome)
	filt(modernGenome)
	ancGenes = utils.myGenomes.Genome(ancGenes, filterOut=filterOut)

#filterOut = set()
#def filt(genome):
#	for (chrom,l) in genome.lstGenes.iteritems():
#		if len(l) < arguments["minimalLength"]:
#			for gene in l:
#				filterOut.update(gene.names)
#filt(ancGenome)
#filt(modernGenome)
#ancGenes = utils.myGenomes.Genome(ancGenes, filterOut=filterOut)

# Renvoie la liste des adjacences entre blocs et des blocs aux extremities
def getAllAdj(genome, dic):

	adj = set()
	extr1 = set()
	extr2 = set()

	for (chrom,l) in genome.lstGenes.iteritems():
		tmp = [(dic[(chrom,i)],gene.strand) for (i,gene) in enumerate(genome.lstGenes[chrom]) if (chrom,i) in dic]
		# Suite de (bloc,orientation)
		tmp = [[(b,s*sa) for (b,sa) in x] for (x,s) in tmp]

		if len(tmp) > 0:

			# Toutes les adjacences
			for (l1,l2) in utils.myTools.myIterator.slidingTuple(tmp):
				for ((b1,s1),(b2,s2)) in itertools.product(l1, l2):
					if (b1,s1) != (b2,s2):
						assert b1 != b2
						adj.add( ((b1,s1),(b2,s2)) )
						adj.add( ((b2,-s2),(b1,-s1)) )

			# Extremites (debut de chrom)
			for (i0,s0) in tmp[0]:
				extr1.add( (i0,s0) )
				extr2.add( (i0,-s0) )

			# Extremites (fin de chrom)
			for (i1,s1) in tmp[-1]:
				extr2.add( (i1,s1) )
				extr1.add( (i1,-s1) )

	return (adj, extr1, extr2)

dicA = collections.defaultdict(list)
dicM = collections.defaultdict(list)
stats = []
dicContent = []

# Parcours des blocs de syntenie
for (n,((ca,da),(cm,dm),_)) in enumerate(myDiags.calcDiags(ancGenome, modernGenome, ancGenes, arguments["minimalLength"]-1, True, myDiags.OrthosFilterType[2])):
	dicContent.append( ((cm,dm[0],dm[-1]), (ca,da[0],da[-1])) )

	# Les scaffolds ne sont pas utilises
	if cm not in modernGenome.chrSet[utils.myGenomes.ContigType.Chromosome]:
		continue
	if len(da) < arguments["minimalLength"]:
		continue
	stats.append(len(da))

	# Dictionnaires pour association gene -> bloc
	for ((im,sm),(ia,sa)) in itertools.izip(dm, da):
		#assert (c1,i1) not in dicM
		dicM[(cm,im)].append( (n,sa) )
		assert (ca,ia) not in dicA
		dicA[(ca,ia)].append( (n,sa) )

for x in dicA.itervalues():
	assert len(x) == 1

print >> sys.stderr, "Synteny:", utils.myMaths.myStats.txtSummary(stats)
blockAdjA = getAllAdj(ancGenome, dicA)
print >> sys.stderr, "anc", len(dicA), len(blockAdjA[0]), len(blockAdjA[1])

# A cause des trous dans les diagonales, les blocs peuvent etre entremeles
# il ne faut garder que les (a,1)-(a+1,1) et (a,-1)-(a-1,-1)
blockAdjA[0].difference_update(set(x for x in blockAdjA[0] if (x[0][1] != x[1][1]) or (x[1][0] != (x[0][0] + x[0][1]))))

blockAdjM = getAllAdj(modernGenome, dicM)
print >> sys.stderr, "mod", len(dicM), len(blockAdjM[0]), len(blockAdjM[1])

def doCmp(blockAdj1, blockAdj2, txt, ind, genome):

	# Le nom du gene a l'extremite
	def getExtr(x, t):
		c = dicContent[x[0]][ind][0]
		(i,_) = dicContent[x[0]][ind][2 if t*x[1] > 0 else 1]
		return genome.lstGenes[c][i].names[0]
	
	# La liste des blocs
	for x in xrange(len(dicContent)):
		print "%s-BLOCK" % txt, x, getExtr((x,1), -1), getExtr((x,1), 1)

	n1 = 0
	n2 = 0
	n3 = 0
	tmp = []
	# Test de chaque adjacence
	for (x1,x2) in blockAdj1[0]:
		if x1 > x2:
			continue
		if (x1,x2) in blockAdj2[0]:
			print "%s-OK" % txt, getExtr(x1, 1), getExtr(x2, -1), x1, x2
			n1 += 1
		elif (x1 in blockAdj2[1]) and (x2 in blockAdj2[2]):
			print "%s-EXTR" % txt, getExtr(x1, 1), getExtr(x2, -1), x1, x2
			n2 += 1
		else:
			print "%s-NO" % txt, getExtr(x1, 1), getExtr(x2, -1), x1, x2
			n3 += 1

	print >> sys.stderr, n1, "OK,", n2, "extremites,", n3, "NO"

print >> sys.stderr, "adjacences de l'ancetre chez le moderne:",
doCmp(blockAdjA, blockAdjM, "anc", 1, ancGenome)

print >> sys.stderr, "adjacences du moderne chez l'ancetre:",
doCmp(blockAdjM, blockAdjA, "mod", 0, modernGenome)

