#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Charge deux fichiers d'arbres et renvoie les differences
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("tree1",file), ("tree2",file)], [], __doc__)

class bufferedTreeReader:
	def __init__(self, filename):
		self.reader = utils.myProteinTree.loadTree(filename)
		self.buffer = []
		self.finished = False

	def readNext(self):
		try:
			self.buffer.append(self.reader.next())
		except StopIteration:
			self.finished = True


reader1 = bufferedTreeReader(arguments["tree1"])
reader2 = bufferedTreeReader(arguments["tree2"])

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



turn = 1
while not (reader1.finished or reader2.finished):
	for ((root1, data1, info1), (root2, data2, info2)) in itertools.product(reader1.buffer, reader2.buffer):
		if root1 == root2:
			lookup(root1)
			reader1.buffer.remove( (root1, data1, info1) )
			reader2.buffer.remove( (root2, data2, info2) )
			break
	else:
		reader1.readNext()
		reader2.readNext()

if len(reader1.buffer) > 0:
	print "roots +1-2", [root1 for (root1,_,_) in reader1.buffer]
	print "roots -1+2", [root2 for (root2,_,_) in reader2.buffer]


