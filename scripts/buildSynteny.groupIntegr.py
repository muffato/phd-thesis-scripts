#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Met bout a bout les diagonales integrees pour former des chromosomes
"""

import sys
import itertools
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths

import utils.myDiags

# Arguments
modesOrthos = list(utils.myDiags.OrthosFilterType._keys)
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("minimalWeight",int,1), ("minimalLength",int,2), ("orthosFilter",str,modesOrthos), \
	("IN.ancGenomeFromDiags",str,""), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


orthosFilter = utils.myDiags.OrthosFilterType[modesOrthos.index(arguments["orthosFilter"])]

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les especes a utiliser
tmp = arguments["target"].split(',')
assert len(arguments["target"]) > 0, "Aucune cible indiquee pour l'extraction des diagonales"

listSpecies = []
for x in tmp:
	listSpecies.extend(phylTree.species[x])

listAncestr = []
for x in tmp:
	listAncestr.extend(set(phylTree.allDescendants[x]).difference(phylTree.species[x]))

dicGenomes = {}
for e in listSpecies:
	dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
genesAnc = {}
for anc in listAncestr:
	dicGenomes[anc] = utils.myGenomes.Genome(arguments["IN.ancGenomeFromDiags"] % phylTree.fileName[anc])
	genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])

def rewriteGenome(genome, dic):
	newGenome = {}
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		tmp = [(dic[(chrom,i)],gene.strand) for (i,gene) in enumerate(genome.lstGenes[chrom]) if (chrom,i) in dic]
		tmp = [(i,s*sa) for ((i,sa),s) in tmp]
		if len(tmp) > 0:
			newGenome[chrom] = [i for (i,l) in itertools.groupby(tmp)]
	return newGenome

def getExtremities(genome):
	extr1 = {}
	extr2 = {}
	for (chrom,l) in genome.iteritems():
		(i0,s0) = l[0]
		(i1,s1) = l[-1]
		extr1[(i0,s0)] = (chrom,1)
		extr2[(i1,s1)] = (chrom,1)
		extr1[(i1,-s1)] = (chrom,-1)
		extr2[(i0,-s0)] = (chrom,-1)
	return (extr1, extr2)

def getAllAdj(anc):
	allAdj = {}
	for esp in listSpecies:

		dicA = {}
		dicM = {}
		stats = []
		print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (anc,esp),
		
		for (n,((c1,d1),(c2,d2),da)) in enumerate(utils.myDiags.calcDiags(dicGenomes[esp], dicGenomes[anc], genesAnc[phylTree.dicParents[anc][esp]], 0, True, orthosFilter)):
			if len(da) < arguments["minimalLength"]:
				continue
			for ((i1,s1),(i2,s2)) in itertools.izip(d1, d2):
				dicM[(c1,i1)] = (n,s1)
				dicA[(c2,i2)] = (n,s1)
			stats.append(len(da))

		print >> sys.stderr, utils.myMaths.myStats.txtSummary(stats),

		newGA = rewriteGenome(dicGenomes[anc], dicA)
		newGM = rewriteGenome(dicGenomes[esp], dicM)
		
		(extr1,extr2) = getExtremities(newGA)
		
		tmp = []
		for (cM,l) in newGM.iteritems():
			for (x1,x2) in utils.myTools.myIterator.slidingTuple(l):
				if (x1 in extr2) and (x2 in extr1):
					tmp.append((extr2[x1],extr1[x2]))

		allAdj[esp] = set(tmp)
		allAdj[esp].update( ((i2,-s2),(i1,-s1)) for ((i1,s1),(i2,s2)) in tmp )
		print >> sys.stderr, "%d adjacences / %d blocs" % (len(tmp), len(newGA))
	return allAdj

toStudy = collections.defaultdict(list)
for (e1,e2) in itertools.combinations(listSpecies, 2):
	for anc in phylTree.dicLinks[e1][e2][1:-1]:
		toStudy[anc].append((e1,e2))

for anc in listAncestr:

	allAdj = getAllAdj(anc)

	gr = utils.myDiags.WeightedDiagGraph()
	for (e1,e2) in toStudy[anc]:
		for x in allAdj[e1] & allAdj[e2]:
			gr.addDiag(x)
	gr.printIniGraph()
	gr.cleanGraphTopDown(2*arguments["minimalWeight"])

	stats = []
	f = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % anc)
	others = set(dicGenomes[anc].lstGenes)
	for (n,(res,scores)) in enumerate(gr.getBestDiags()):
		i = 0
		for (x,s) in res:
			others.remove(x)
			if s > 0:
				for gene in dicGenomes[anc].lstGenes[x]:
					f.csvobject.writerow( [n, i, i, gene.strand, " ".join(gene.names)] )
					i += 1
			else:
				for gene in dicGenomes[anc].lstGenes[x].__reversed__():
					f.csvobject.writerow( [n, i, i, -gene.strand, " ".join(gene.names)] )
					i += 1
		stats.append(i)
	
	for x in others:
		for (i,gene) in enumerate(dicGenomes[anc].lstGenes[x]):
			f.csvobject.writerow( [n, i, i, gene.strand, " ".join(gene.names)] )
		n += 1
		if i > 0:
			stats.append(i+1)

	print >> sys.stderr, anc, utils.myMaths.myStats.txtSummary(stats)
	f.file.close()

