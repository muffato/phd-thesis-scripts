#!/usr/bin/env python2

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

	f = utils.myFile.openFile(arguments["diags"] % phylTree.fileName[anc], "r")
	for l in f:
		t = l.split("\t")
		lg = [int(x) for x in t[2].split()]
		ls = [int(x) for x in t[3].split()]
		for ((g1,s1),(g2,s2)) in utils.myTools.myIterator.slidingTuple(zip(lg,ls)):
			pos1 = [phylTree.dicGenes[x] for x in ancGenes.lstGenes[None][g1].names if x in phylTree.dicGenes]
			pos2 = [phylTree.dicGenes[x] for x in ancGenes.lstGenes[None][g2].names if x in phylTree.dicGenes]
			conserved = set()
			scaff = set()
			notseen = set(phylTree.species[anc])
			for ((e1,c1,i1),(e2,c2,i2)) in itertools.product(pos1, pos2):
				if e1 != e2:
					continue
				notseen.discard(e1)
				if (c1 != c2) and ((i1 == 0) or (i1 == len(phylTree.dicGenomes[e1].lstGenes[c1])-1)) and ((i2 == 0) or (i2 == len(phylTree.dicGenomes[e2].lstGenes[c2])-1)):
					scaff.add(e1)
					continue
				if c1 != c2:
					continue
				if (i2 == i1+1) and (phylTree.dicGenomes[e1].lstGenes[c1][i1].strand == s1) and (phylTree.dicGenomes[e2].lstGenes[c2][i2].strand == s2):
					conserved.update(phylTree.dicLinks[e1][anc])
				elif (i2 == i1-1) and (phylTree.dicGenomes[e1].lstGenes[c1][i1].strand == -s1) and (phylTree.dicGenomes[e2].lstGenes[c2][i2].strand == -s2):
					conserved.update(phylTree.dicLinks[e2][anc])
			lst = [e for e in phylTree.allDescendants[anc] if (e not in conserved) and (phylTree.parent[e][0] in conserved) and (len(phylTree.species[e].difference(scaff)) > 0) and (len(phylTree.species[e].difference(notseen)) > 0)]
			print utils.myFile.myTSV.printLine([anc, g1, g2, ancGenes.lstGenes[None][g1].names[0], ancGenes.lstGenes[None][g2].names[0], s1, s2, len(lst), "|".join(lst)])
			#print utils.myFile.myTSV.printLine([anc, g1, g2, ancGenes.lstGenes[None][g1].names[0], ancGenes.lstGenes[None][g2].names[0], s1, s2, len(lst), "|".join(lst), "|".join(e for e in phylTree.allDescendants[anc] if (e not in conserved)), "|".join(e for e in phylTree.allDescendants[anc] if (e not in conserved) and (phylTree.parent[e][0] in conserved)), "|".join(conserved), "|".join(scaff)])
	f.close()

