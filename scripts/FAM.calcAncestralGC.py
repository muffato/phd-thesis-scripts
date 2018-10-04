#!/usr/bin/env python2

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [], [("range",str,""), ("phylTree",str,""), ("alignment-FASTA",str,"")], "Retrouve le GC ancestral" )

for treeID in utils.myTools.getRange(arguments["range"]):
	
	print >> sys.stderr, treeID, "...",

	tree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"] % treeID)

	# Chargement des CDS
	seq = utils.myGenomes.myFASTA.loadFile(arguments["alignment-FASTA"] % treeID)

	res = []
	n = len(seq.values()[0])
	for i in xrange(n):
		values = {}
		for (e,s) in seq.iteritems():
			if i < len(s):
				c = s[i]
			else:
				print >> sys.stderr, "!1", e, s, i, n
				c = "-"
			if c in "GCgc":
				values[e] = 1
			elif c in "ATat":
				values[e] = 0
		res.append(tree.calcWeightedValue(values, -1, None))
	for (ie,e) in enumerate(tree.allNames):
		l = [r[ie] for r in res if r[ie] >= 0]
		print utils.myFile.myTSV.printLine([treeID, " ".join(tree.species[e]), len(l), utils.myMaths.myStats.myStats.mean(l)])

