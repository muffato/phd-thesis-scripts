#! /users/ldog/muffato/python

__doc__ = """
	Parcourt les deux genomes et extrait les intervalles conserves, les points de cassure
	  en tenant compte des evenements de genes (gain/perte/duplication)
"""

import sys
import collections
import itertools

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( [("phylTree.conf",file)], [("modernGenomes",str,""), ("ancGenomes",str,""), ("ancGenesFile",str,""), ("diags",str,"")], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

def rewriteGenome(genome, dic):
	newGenome = {}
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		tmp = [(dic[gene.names[0]],gene.strand) for gene in genome.lstGenes[chrom] if gene.names[0] in dic]
		tmp = [(i,s*sa) for ((i,sa),s) in tmp]
		if len(tmp) > 0:
			newGenome[chrom] = [i for (i,l) in itertools.groupby(tmp)]
	return newGenome

def getAdj(genome):
	adj = set()
	extr1 = set()
	extr2 = set()
	for (chrom,l) in genome.iteritems():
		adj.update(utils.myTools.myIterator.slidingTuple(l))
		extr1.add( l[0] )
		extr2.add( l[-1] )
		extr1.add( (l[-1][0],-l[-1][1]) )
		extr2.add( (l[0][0],-l[0][1]) )
	return (adj,extr1,extr2)

def check(adjA, adjM, extrM1, extrM2, txt, dic):
	def span(l):
		return "%s >%d> %s" % (l[0], len(l), l[-1])
	for ((i,si),(j,sj)) in adjA:
		if (((i,si),(j,sj)) not in adjM) and (((j,-sj),(i,-si)) not in adjM):
			if ((i,si) not in extrM2) or ((j,sj) not in extrM1):
				#print txt, (i,j), dic[i], dic[j]
				print txt, "%d/%d" % (i,si), "+", "%d/%d" % (j,sj), ":", span(dic[i][0]), "/", span(dic[i][1]), "!=!", span(dic[j][0]), "/", span(dic[j][1])


allgenomes = {}
def load(esp):
	if esp not in allgenomes:
		if esp in phylTree.listSpecies:
			allgenomes[esp] = utils.myGenomes.Genome(arguments["modernGenomes"] % phylTree.fileName[esp])
		else:
			allgenomes[esp] = utils.myGenomes.Genome(arguments["ancGenomes"] % phylTree.fileName[esp], ancGenes=utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[esp]))
	return allgenomes[esp]


# Chargement des fichiers de diagonales valides
for anc in phylTree.listAncestr:
	if not utils.myFile.hasAccess(arguments["diags"] % phylTree.fileName[anc]):
		continue

	dicA = collections.defaultdict(dict)
	dicM = collections.defaultdict(dict)
	dicD = collections.defaultdict(dict)
	f = utils.myFile.openFile(arguments["diags"] % phylTree.fileName[anc], "r")
	for (i,l) in enumerate(f):
		t = l[:-1].split("\t")
		if phylTree.isChildOf(t[2], t[5]) and (t[5] == anc):
			dS = [int(x) for x in t[9].split()]
			dA = t[7].split()
			dM = t[4].split()
			for (x,s) in itertools.izip(dA, dS):
				dicA[t[2]][x] = (i,s)
			for (x,s) in itertools.izip(dM, dS):
				dicM[t[2]][x] = (i,s)
			dicD[t[2]][i] = (dA, dM)
	f.close()
	if len(dicD) == 0:
		continue

	ancGenome = load(anc)

	def searchDesc(esp):
		if esp in dicD:
			genome = load(esp)

			newGA = rewriteGenome(ancGenome, dicA[esp])
			newGM = rewriteGenome(genome, dicM[esp])

			#for cA in newGA:
			#	print >> sys.stderr, "ANC", cA, newGA[cA]
			#for cM in newGM:
			#	print >> sys.stderr, "MOD", cM, newGM[cM]

			(adjA,extrA1,extrA2) = getAdj(newGA)
			(adjM,extrM1,extrM2) = getAdj(newGM)

			check(adjA, adjM, extrM1, extrM2, "LOST %s/%s" % (anc,esp), dicD[esp])
			check(adjM, adjA, extrA1, extrA2, "GAIN %s/%s" % (anc,esp), dicD[esp])
		elif esp in phylTree.items:
			for (e,_) in phylTree.items[esp]:
				searchDesc(e)

	searchDesc(anc)
