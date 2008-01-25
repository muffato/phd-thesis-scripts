#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( [], [("start",int,0), ("end",int,0), ("phylTree",str,""), ("alignment-FASTA",str,"")], "Retrouve le GC ancestral" )

for treeID in xrange(options["start"], options["end"]+1):
	
	print >> sys.stderr, treeID, "...",

	tree = utils.myPhylTree.PhylogeneticTree(options["phylTree"] % treeID, buildLinks = False)

	# Chargement des CDS
	seq = utils.myGenomes.loadFastaFile(options["alignment-FASTA"] % treeID)

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
		res.append(tree.calcWeightedValue(values, -1, None, None))
	for (ie,e) in enumerate(tree.allNames):
		l = [r[ie] for r in res if r[ie] >= 0]
		print utils.myTools.printLine([treeID, " ".join(tree.species[e]), len(l), utils.myMaths.mean(l)])

