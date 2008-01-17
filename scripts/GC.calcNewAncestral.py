#! /users/ldog/muffato/python -OO

import sys
import utils.myMaths
import utils.myTools
import utils.myPhylTree

(noms_fichiers, options) = utils.myTools.checkArgs( [], [("start",int,0), ("end",int,0), ("phylTree",str,""), ("gc3.mase",str,"")], "Retrouve le GC ancestral" )

for treeID in xrange(options["start"], options["end"]+1):
	
	print >> sys.stderr, treeID

	tree = utils.myPhylTree.PhylogeneticTree(options["phylTree"] % treeID, buildLinks = False)

	f = utils.myTools.myOpenFile(options["gc3.mase"] % treeID, "r")
	nom = None
	seq = {}
	for l in f:
		if l[0] == ";":
			continue
		if nom == None:
			nom = l[:-1]
		else:
			seq[nom] = l[:-1]
			n = len(seq[nom])
			nom = None
	f.close()

	res = []
	for i in xrange(n):
		values = {}
		for s in seq:
			if i < len(seq[s]):
				c = seq[s][i]
			else:
				print >> sys.stderr, "!1", s, i, n
				c = "-"
			if c in "GCgc":
				values[s] = 1
			elif c in "ATat":
				values[s] = 0
		res.append(tree.calcWeightedValue(values, -1, None, None))
	for (ie,e) in enumerate(tree.allNames):
		l = [r[ie] for r in res if r[ie] >= 0]
		print utils.myTools.printLine([treeID, tree.species[e], len(l), utils.myMaths.mean(l)])

