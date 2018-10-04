#!/usr/bin/env python2

__doc__ = """
	Cree la base Genomicus (table blocks_stats)
	Necessite la table Gene pour positionner les blocs sur les chromosomes
	Lit les fichiers de diagonales + genes ancestraux
"""

import sys
import collections

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree



def storeData(genome, esp):

	contigs = [len(x) for x in genome.lstGenes.itervalues() if len(x) >= 2]
	print utils.myFile.myTSV.MySQLFileWriter( (
		phylTree.indNames[esp],
		sum(len(x) for x in genome.lstGenes.itervalues()),
		len(contigs),
		sum(contigs),
		utils.myMaths.myStats.getValueNX(sorted(len(x) for x in genome.lstGenes.itervalues()), 50)
		))
	print >> sys.stderr, "OK"


arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("modernGenomes",str,""), ("ancGenes",str,""), ("diags",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for esp in phylTree.listSpecies:
	genome = utils.myGenomes.Genome(arguments["modernGenomes"] % phylTree.fileName[esp])
	storeData(genome, esp)
	

for anc in phylTree.listAncestr:
	ancGenes = utils.myGenomes.Genome(arguments["ancGenes"] % phylTree.fileName[anc])
	genome = utils.myGenomes.Genome(arguments["diags"] % phylTree.fileName[anc], ancGenes=ancGenes)
	storeData(genome, anc)



