#!/usr/bin/env python2

__doc__ = """
	Renvoie le nombre d'orthologues vs le nombre d'exons
"""

import sys
import collections

import utils.myFile
import utils.myGenomes
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("nbExonsFile",file), ("ancestr",str)],\
	[("genesFile",str,""), ("ancGenesFile",str,"")], \
	__doc__ \
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
phylTree.loadAllSpeciesSince(arguments["ancestr"], arguments["genesFile"], storeGenomes=False)

ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % arguments["ancestr"])

nbExons = {}
f = utils.myFile.openFile(arguments["nbExonsFile"], "r")
for l in f:
	t = l.split()
	nbExons[t[0]] = int(t[1])
f.close()

for g in ancGenes:
	dicEsp = collections.defaultdict(list)
	for x in g.names:
		if x in phylTree.dicGenes:
			if phylTree.dicGenes[x].species not in phylTree.species["Euarchontoglires"]:
				continue
			dicEsp[phylTree.dicGenes[x].species].append(x)
	lgh = dicEsp["Homo sapiens"]
	if len(lgh) == 0:
		continue
	lgh1 = [x for x in lgh if nbExons[x] == 1]
	lghN = [x for x in lgh if nbExons[x] > 1]

	le = [e for e in dicEsp if e not in phylTree.lstEsp2X]

	print utils.myFile.myTSV.printLine([len(lgh1), " ".join(lgh1), len(lghN), " ".join(lghN), len(dicEsp)-1, sum(len(dicEsp[x]) for x in dicEsp)-len(lgh), len(le)-1, sum(len(dicEsp[x]) for x in le)-len(lgh)])


