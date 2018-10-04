#!/usr/bin/env python2

__doc__ = """
	Cherche des points de rearrangements par syntenie homme/ancetre
"""

import sys
import operator
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

# Chargement des genomes
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

dicM = {}
stats = []
dicContent = []
genomeA = collections.defaultdict(list)

# Parcours des blocs de syntenie (selon le genome ancestral)
for (n,((ca,da),(cm,dm),_)) in enumerate(myDiags.calcDiags(ancGenome, modernGenome, ancGenes, 0, True, myDiags.OrthosFilterType[2])):
#for (n,((ca,da),(cm,dm),_)) in enumerate(myDiags.calcDiags(ancGenome, modernGenome, ancGenes, arguments["minimalLength"]-1, True, myDiags.OrthosFilterType[2])):
	
	dicContent.append( ( \
		(cm, modernGenome.lstGenes[cm][dm[0][0]].names[0], modernGenome.lstGenes[cm][dm[-1][0]].names[0]), \
		(ca, ancGenome.lstGenes[ca][da[0][0]].names[0], ancGenome.lstGenes[ca][da[-1][0]].names[0]) \
	) )
	print "BLOCK %d:" % n, len(da),
	print "|", utils.myFile.myTSV.printLine(dicContent[-1][1], delim=" "),
	print "|", utils.myFile.myTSV.printLine(dicContent[-1][0], delim=" ")

	# Les scaffolds ne sont pas utilises
	if cm not in modernGenome.chrSet[utils.myGenomes.ContigType.Chromosome]:
		continue
	
	# Les blocs trop courts sont oublies
	if len(da) < arguments["minimalLength"]:
		continue
	stats.append(len(da))

	genomeA[ca].append(n)
	
	# Dictionnaires pour association gene -> bloc
	for ((im,_),(_,sa)) in itertools.izip(dm, da):
		assert (cm,im) not in dicM
		dicM[(cm,im)] = (n,sa)

print >> sys.stderr, "synteny:", utils.myMaths.myStats.txtSummary(stats), len(dicM), "genes"

AdjStruct = collections.namedtuple("AdjStruct", ["adj", "extr1", "extr2"])

# Met a jour une structure d'adjacences, a partir du contenu en blocs d'un chromosome
def updateAdj(allAdj, chrom):
	
	# Toutes les adjacences
	for ((b1,s1),(b2,s2)) in utils.myTools.myIterator.slidingTuple(chrom):
		allAdj.adj.add( ((b1,s1),(b2,s2)) )
		allAdj.adj.add( ((b2,-s2),(b1,-s1)) )
	
	# Extremites (debut de chrom)
	(i0,s0) = chrom[0]
	allAdj.extr1.add( (i0,s0) )
	allAdj.extr2.add( (i0,-s0) )

	# Extremites (fin de chrom)
	(i1,s1) = chrom[-1]
	allAdj.extr2.add( (i1,s1) )
	allAdj.extr1.add( (i1,-s1) )


# Construction des adjacences de l'ancetre
allAdjA = AdjStruct(set(), set(), set())
for chrom in genomeA:
	updateAdj(allAdjA, [(i,1) for i in genomeA[chrom]])
print >> sys.stderr, "anc:", len(allAdjA.adj)/2, "adjacences,", len(allAdjA.extr1), "extremites"

# Construction des adjacences de l'espece moderne
allAdjM = AdjStruct(set(), set(), set())
for (chrom,l) in modernGenome.lstGenes.iteritems():

	# Les genes presents dans des blocs
	tmp = [(dicM[(chrom,i)],gene.strand) for (i,gene) in enumerate(l) if (chrom,i) in dicM]
	# Liste de (bloc,orientation)
	tmp = [(b,s*sa) for ((b,sa),s) in tmp]

	if len(tmp) == 0:
		continue

	# Les blocs peuvent se chevaucher mais pas s'inclure
	# On peut donc les ordonner
	orderedBlocks = []
	seenBlocks = set()
	for (i,b) in enumerate(tmp):
		if b not in seenBlocks:
			seenBlocks.add(b)
			orderedBlocks.append(b)

	updateAdj(allAdjM, orderedBlocks)

print >> sys.stderr, "mod:", len(allAdjM.adj)/2, "adjacences,", len(allAdjM.extr1), "extremites"

def doCmp(allAdjT, allAdjR, txt):

	# Extremites des blocs, en fonction des orientations
	def printExtr((i1,s1), (i2,s2)):
		print dicContent[i1][1][2 if s1 > 0 else 1], dicContent[i2][1][1 if s2 > 0 else 1],
		print dicContent[i1][0][2 if s1 > 0 else 1], dicContent[i2][0][1 if s2 > 0 else 1],
		print "%d/%d" % (i1,s1), "%d/%d" % (i2,s2)

	n1 = 0
	n2 = 0
	n3 = 0
	# Test de chaque adjacence chez l'espece reference
	for (x1,x2) in allAdjT.adj:
		if x1 > x2:
			continue
		if (x1,x2) in allAdjR.adj:
			# Adjacence presente
			print txt + "-ADJ=OK",
			n1 += 1
		elif (x1 in allAdjR.extr1) and (x2 in allAdjR.extr2):
			# Fusion de telomeres
			print txt + "-ADJ=EXTR",
			n2 += 1
		else:
			# Adjacence differente
			print txt + "-ADJ=NO",
			n3 += 1
		printExtr(x1, x2)

	print >> sys.stderr, n1+n2+n3, "adjacences", "=", n1, "OK", "+", n2, "extremites", "+", n3, "NO"

print >> sys.stderr, "anc -> mod:",
doCmp(allAdjA, allAdjM, "ANC")

print >> sys.stderr, "mod -> anc:",
doCmp(allAdjM, allAdjA, "MOD")

