#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Prend deux listes d'intervalles et extrait les intervalles conserves 
"""

import sys
import math
import itertools
import collections

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( [("phylTree",file)], [("ancGenes",str,""), ("genes",str,""), ("diags",str,"")], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])
phylTree.loadAllSpeciesSince(None, arguments["genes"])

for anc in phylTree.listAncestr:
	ancGenes = utils.myGenomes.Genome(arguments["ancGenes"] % phylTree.fileName[anc])

	conserved = collections.defaultdict(set)
	
	#count = {}
	#count[-3] = count[3] = collections.defaultdict(int)
	#count[1] = collections.defaultdict(int)
	#count[-1] = collections.defaultdict(int)

	f = utils.myFile.openFile(arguments["diags"] % phylTree.fileName[anc], "r")
	for l in f:
		t = l.split("\t")
		lg = [int(x) for x in t[2].split()]
		ls = [int(x) for x in t[3].split()]
		for ((g1,s1),(g2,s2)) in utils.myTools.myIterator.slidingTuple(zip(lg,ls)):
			pos1 = [phylTree.dicGenes[x] for x in ancGenes.lstGenes[None][g1].names if x in phylTree.dicGenes]
			pos2 = [phylTree.dicGenes[x] for x in ancGenes.lstGenes[None][g2].names if x in phylTree.dicGenes]
			for ((e1,c1,i1),(e2,c2,i2)) in itertools.product(pos1, pos2):
				if e1 != e2:
					continue
				if c1 != c2:
					continue
				if (i2 == i1+1) and (phylTree.dicGenomes[e1].lstGenes[c1][i1].strand == s1) and (phylTree.dicGenomes[e2].lstGenes[c2][i2].strand == s2):
					conserved[e1].add((c1,i1,i2))
					#conserved[e1].append(phylTree.dicGenomes[e2].lstGenes[c2][i2].beginning-phylTree.dicGenomes[e1].lstGenes[c1][i1].end+1)
					#count[2*s1+s2][e1] += 1
				elif (i2 == i1-1) and (phylTree.dicGenomes[e1].lstGenes[c1][i1].strand == -s1) and (phylTree.dicGenomes[e2].lstGenes[c2][i2].strand == -s2):
					conserved[e1].add((c1,i2,i1))
					#conserved[e1].append(phylTree.dicGenomes[e1].lstGenes[c1][i1].beginning-phylTree.dicGenomes[e2].lstGenes[c2][i2].end+1)
					#count[2*s1+s2][e1] += 1
	f.close()

	for esp in phylTree.species[anc]:
		for c in phylTree.dicGenomes[esp].chrList[utils.myGenomes.ContigType.Chromosome] + phylTree.dicGenomes[esp].chrList[utils.myGenomes.ContigType.Scaffold]:
			chrom = phylTree.dicGenomes[esp].lstGenes[c]
			for i in xrange(len(chrom)-1):
				print utils.myFile.myTSV.printLine([ anc, esp, chrom[i], chrom[i+1], chrom[i+1].beginning-chrom[i].end+1, (c,i,i+1) in conserved[esp] ])
				#print anc, esp, chrom[i].strand, chrom[i+1].strand, chrom[i+1].beginning-chrom[i].end+1, (c,i,i+1) in conserved[esp]
		
		#print utils.myFile.myTSV.printLine([ anc, esp, count[3][esp]+count[1][esp]+count[-1][esp], count[3][esp], count[1][esp], count[-1][esp] ])

