#!/usr/bin/env python2

__doc__ = """
	Compte le nombre de comparaisons faites pour chaque ancetre, et le correle a la longueur des diagonales
"""


import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], [("diags",str,"")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for anc in phylTree.listAncestr:
	n0 = len(phylTree.outgroupSpecies[anc])
	l = [len(phylTree.species[x]) for (x,_) in phylTree.items[anc]]
	l.append(n0)
	r = []
	f = utils.myFile.openFile(arguments["diags"] % phylTree.fileName[anc], "r")
	for line in f:
		x = int(line.split()[1])
		if x > 1:
			r.append(x)
	f.close()
	nbc = sum(n1*n2 for (n1,n2) in utils.myTools.myIterator.tupleOnStrictUpperList(l))
	lll = float(sum(r))/len(r)
	print utils.myFile.myTSV.printLine([anc, nbc, phylTree.ages[anc], lll, float(nbc)/phylTree.ages[anc]])

