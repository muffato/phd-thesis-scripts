#!/usr/bin/env python2

__doc__ = """
	A partir de diagonales pair-wise et de blocs ancestraux fixes,
	   construit des versions integrees qui correspondent a des segments de chromosomes ancestraux
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
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.ancDiags",str,""), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

dicGenomes = {}
for esp in phylTree.listSpecies:
	dicGenomes[esp] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[esp])
ancGenes = {}
for anc in phylTree.listAncestr:
	ancGenes[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	dicGenomes[anc] = utils.myGenomes.Genome(arguments["IN.ancDiags"] % phylTree.fileName[anc], ancGenes=ancGenes[anc])

def do(node):

	def genetype(x):
		if x is None:
			return "N"
		n = count[x]
		assert n >= 1
		if n == 1:
			return "S"
		else:
			return "D %d (%s/%d)" % (n,x[0],x[1])

	genome1 = dicGenomes[node]
	for (f,_) in phylTree.items.get(node, []):
		genome2 = dicGenomes[f]
		newGenome = {}
		count = collections.defaultdict(int)
		for chrom in genome2.lstGenes:
			newGenome[chrom] = [genome1.dicGenes.get(gene.names[-1]) for gene in genome2.lstGenes[chrom]]
			for x in newGenome[chrom]:
				if x is not None:
					count[x] += 1
		for chrom in genome2.lstGenes:
			tmp = [genetype(x) for x in newGenome[chrom]]
			print tmp
			for (t,l) in itertools.groupby(tmp):
				if t != "S":
					print utils.myFile.myTSV.printLine( (f, t.split()[0], len(list(l))) )

		# Dans l'autre sens: pertes
		for chrom in genome1.lstGenes:
			tmp = [genome2.dicGenes.get(gene.names[-1]) for gene in genome1.lstGenes[chrom]]
			for (t,l) in itertools.groupby(tmp):
				if t is None:
					print utils.myFile.myTSV.printLine( (f, "X", len(list(l))) )
		do(f)

do(phylTree.root)

