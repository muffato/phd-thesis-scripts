#!/usr/bin/env python2

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
arguments = utils.myTools.checkArgs( [("ancGenome",file), ("modernGenome",file)], [], __doc__)

ancGenome = utils.myGenomes.Genome(arguments["ancGenome"])
genome = utils.myGenomes.Genome(arguments["modernGenome"])

for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
	
	dataGenes = []
	for gene in genome.lstGenes[chrom]:
		if gene.names[0] in ancGenome.dicGenes:
			(c,_) = ancGenome.dicGenes[gene.names[0]]
			if len(ancGenome.lstGenes[c]) == 1:
				dataGenes.append('S')
			else:
				dataGenes.append('B')
		else:
			dataGenes.append('N')
	
	dataInterv = []
	for (g1,g2) in utils.myTools.myIterator.slidingTuple(genome.lstGenes[chrom]):
		if (g1.names[0] in ancGenome.dicGenes) and (g2.names[0] in ancGenome.dicGenes):
			(c1,i1) = ancGenome.dicGenes[g1.names[0]]
			s1 = ancGenome.lstGenes[c1][i1].strand
			(c2,i2) = ancGenome.dicGenes[g2.names[0]]
			s2 = ancGenome.lstGenes[c2][i2].strand
			if c1 == c2:
				if (i1 == i2) and (g1.strand == g2.strand):
					dataInterv.append('D')
				elif (i2-i1 == g1.strand*s1) and (g1.strand*g2.strand == s1*s2):
					dataInterv.append('S')
				else:
					dataInterv.append('R')
			else:
				extr1 = ((g1.strand == s1) and (i1 == len(ancGenome.lstGenes[c1])-1)) or ((g1.strand == -s1) and (i1 == 0))
				extr2 = ((g2.strand == -s2) and (i2 == len(ancGenome.lstGenes[c2])-1)) or ((g2.strand == s2) and (i2 == 0))
				if extr1 and extr2:
					dataInterv.append('?')
				else:
					dataInterv.append('R')
		else:
			dataInterv.append('N')
	dataInterv = ['*'] + [x.lower() for x in dataInterv] + ['*']

	for i in xrange(len(genome.lstGenes[chrom])-1):
		print chrom, i, genome.lstGenes[chrom][i].names[0], genome.lstGenes[chrom][i+1].names[0], dataInterv[i], dataGenes[i], dataInterv[i+1], dataGenes[i+1], dataInterv[i+2]


