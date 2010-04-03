#! /users/ldog/muffato/python

__doc__ = """
	Charge deux fichiers d'arbres et renvoie les differences
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("tree1",file), ("tree2",file)], [], __doc__)


roots1 = set()
data1 = {}
info1 = {}
for (root,data,info) in utils.myProteinTree.loadTree(arguments["tree1"]):
	roots1.add(root)
	data1.update(data)
	info1.update(info)

roots2 = set()
data2 = {}
info2 = {}
for (root,data,info) in utils.myProteinTree.loadTree(arguments["tree2"]):
	roots2.add(root)
	data2.update(data)
	info2.update(info)

def sameinfo(node, inf1, inf2):
	if inf1 != inf2:
		inf1 = set(inf1.items())
		inf2 = set(inf2.items())
		print "info +1-2", node, inf1.difference(inf2)
		print "info -1+2", node, inf2.difference(inf1)

def lookup(node):
	if (node in data1) ^ (node in data2):
		print "anc vs species:", node, node in data1, node in data2
	elif node in data1:
		sameinfo(node, info1[node], info2[node])
		d1 = set(data1[node])
		d2 = set(data2[node])
		if d1 != d2:
			print "data +1-2", node, d1.difference(d2)
			print "data -1+2", node, d2.difference(d1)
		for (x,_) in (d1 & d2):
			lookup(x)
	else:
		sameinfo(node, info1[node], info2[node])

if roots1 != roots2:
	print "roots +1-2", roots1.difference(roots2)
	print "roots -1+2", roots2.difference(roots1)
for r in (roots1 & roots2):
	lookup(r)

