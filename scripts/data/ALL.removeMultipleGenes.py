#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Analyse l'arbre et supprime les genes qui sont presents via plusieurs transcrits
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("ensemblTree",file)], [], __doc__)

count = collections.defaultdict(list)

allTrees = {}
allRoots = []

for (root,data,info) in utils.myProteinTree.loadTree(arguments["ensemblTree"]):

	allTrees[root] = (data,info)
	allRoots.append(root)

	def countGenes(node):

		if node in data:
			for (g,_) in data[node]:
				countGenes(g)
		else:
			count[(info[node]['gene_name'],info[node]['taxon_name'])].append((node,root))

	countGenes(root)


for ((name,esp),l) in count.iteritems():
	if len(l) == 1:
		continue
	delete = True
	s = set(root for (node,root) in l)
	print >> sys.stderr, name, esp, len(l), len(s),
	if len(s) == 1:
		root = s.pop()
		lca = None
		(data,info) = allTrees[root]
		def findLastCommonAncestor(node):
			global lca
			global delete
			if node in data:
				nOK = nTot = 0
				for (g,_) in data[node]:
					res = findLastCommonAncestor(g)
					nOK += res[0]
					nTot += res[1]
				if (nOK == nTot) and (nOK == len(l)):
					lca = node
					delete = False
				return (nOK, nTot)

			else:
				if info[node]['gene_name'] == name:
					return (1,1)
				else:
					return (0,1)
		findLastCommonAncestor(root)
	
	print >> sys.stderr, delete
	# Rebuild the tree
	if delete:
		for (node,root) in l:
			(data,info) = allTrees[root]
			def recdelete(x):
				if x in data:
					ll = [(g,l) for (g,l) in data[x] if recdelete(g)]
					if len(ll) == 0:
						del data[x]
					else:
						data[x] = ll
					return x in data
					#data[x] = [(g,l) for (g,l) in data[x] if g != node]
					#for (g,_) in data[x]:
					#	recdelete(g)
				else:
					return x != node
			recdelete(root)
	else:
		del data[lca]
		info[lca]['taxon_name'] = esp
		info[lca]['Duplication'] = 0
		info[lca]['gene_name'] = name
		info[lca]['transcript_name'] = "/".join(info[node]['transcript_name'] for (node,_) in l)
		info[lca]['protein_name'] = "/".join(info[node]['protein_name'] for (node,_) in l)


for root in allRoots:
	(data,info) = allTrees[root]
	utils.myProteinTree.printTree(sys.stdout, data, info, root)

