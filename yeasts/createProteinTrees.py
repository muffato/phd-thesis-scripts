#!/usr/bin/env python2

__doc__ = """
	Convertit les arbres de REGEV au bon format
"""


import sys
import cStringIO
import collections

import utils.myPhylTree
import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs( [("speciesTree",file), ("all-output.txt",file)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

nodeid = 0


def processData(data):

	tree = utils.myPhylTree.PhylogeneticTree(cStringIO.StringIO(data[2]))
	
	@utils.myTools.memoize
	def lastCommonAncestor(node):
		if node in tree.items:
			return phylTree.lastCommonAncestor([lastCommonAncestor(e) for (e,_) in tree.items[node]])
		else:
			return phylTree.officialName[node.split("|")[0]]

	@utils.myTools.memoize
	def isDupNode(node):
		if node in tree.items:
			return lastCommonAncestor(node) in [lastCommonAncestor(e) for (e,_) in tree.items[node]]
		else:
			return False
	
	def printTree(indent, node):
		global nodeid
		#print "%sid\t%d" % (indent, nodeid)
		nodeid += 1
		count[lastCommonAncestor(node)] += 1
		dup[lastCommonAncestor(node)] += int(isDupNode(node))
		
		info = {"taxon_name": lastCommonAncestor(node), "Duplication": 2*int(isDupNode(node))}
		if node not in tree.items:
			x = node.split("|")
			info["gene_name"] = x[1]
			#info["taxon_name"] = [y for y in phylTree.commonNames[x[0]] if not isinstance(y, int) and not y.startswith("node") and len(y) > 5][0]

		#print "%sinfo\t%s" % (indent, info)

		if node in tree.items:
			indent = indent + "\t"
			for (e,l) in tree.items[node]:
				#print "%slen\t%g" % (indent,l)
				printTree(indent, e)

	count = collections.defaultdict(int)
	dup = collections.defaultdict(int)
	printTree("", tree.root)

	print data[5]
	print data[6]
	print "New Counts:\t " + " ".join(str(count["node"+str(i)]) for i in xrange(45))
	print data[7]
	print "New Dup:\t " + " ".join(str(dup["node"+str(i)]) for i in xrange(45))
	print
	
	#counts = sum(int(x) for x in data[6].split("\t")[1].split())
	#dups = sum(int(x) for x in data[7].split("\t")[1].split())
	#assert dups >= len([x for x in tree.listAncestr if isDupNode(x)]), data
	#assert (counts + dups) >= (len(tree.listSpecies) + len(tree.listAncestr)), (counts, dups, len(tree.listSpecies), len(tree.listAncestr))



currData = None
f = utils.myFile.openFile(arguments["all-output.txt"], "r")
for l in f:
	if l.startswith(">"):
		if currData is not None:
			processData(currData)
		currData = []
	currData.append(l[:-1])
f.close()
if currData is not None:
	processData(currData)

