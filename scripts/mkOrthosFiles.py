#!/usr/bin/env python2

import sys
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("dicGeneProtTrans",file)], \
	[("genesFile",str,""), ("ancGenesFile",str,""), ("outputFile",str,"")], \
	"Cree les fichiers d'orthologues pair-wise" \
)



phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
phylTree.loadAllSpeciesSince(None, arguments["genesFile"], storeGenomes=False)
ancGenes = {}
for anc in phylTree.listAncestr:
	ancGenes[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])

# Association gene_name > transcript_name
f = utils.myFile.openFile(arguments["dicGeneProtTrans"], "r")
dicProtTransc = {}
for l in f:
	t = l.split()
	dicProtTransc[t[0]] = t[2]
f.close()


#for (e1,e2) in itertools.product(["Homo sapiens"], phylTree.listSpecies):
for (e1,e2) in itertools.product(phylTree.listSpecies, phylTree.listSpecies):
	if e1 == e2:
		continue
	anc = phylTree.dicParents[e1][e2]
	print >> sys.stderr, "%s: %s/%s ..." % (anc,e1,e2),
	f = utils.myFile.openFile(arguments["outputFile"] % (phylTree.fileName[e1],phylTree.fileName[e2]), "w")
	for g in ancGenes[anc]:
		l1 = [dicProtTransc[x] for x in g.names[1:] if phylTree.dicGenes[x][0]==e1 and x in dicProtTransc]
		l2 = [dicProtTransc[x] for x in g.names[1:] if phylTree.dicGenes[x][0]==e2 and x in dicProtTransc]
		#l1 = [dic[e1][x] for x in g.names[1:] if x in dic[e1]]
		#l2 = [dic[e2][x] for x in g.names[1:] if x in dic[e2]]
		if (len(l1) > 0) and (len(l2) > 0):
			print >> f, len(l1), len(l2), " ".join(l1), " ".join(l2)
	f.close()
	print >> sys.stderr, "OK"


