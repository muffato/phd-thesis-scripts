#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Prend deux listes d'intervalles et extrait les intervalles conserves 
"""

import sys
import math
import itertools
import collections

import utils.myFile
import utils.myDiags
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Argument:
arguments = utils.myTools.checkArgs( [("list1",file), ("list2",file), ("ancGenes",file)], [], __doc__)

ancGenes = utils.myGenomes.Genome(arguments["ancGenes"])

def loadData(filename):
	data = collections.defaultdict(list)
	back = []
	f = utils.myFile.openFile(filename, "r")
	for l in f:
		l = l[:-1]
		t = l.split("\t")
		ga = tuple(ancGenes.dicGenes[g][1] for g in t[1].split('/') if g in ancGenes.dicGenes)
		if len(frozenset(ga)) == 2:
			gs = tuple(int(s) for s in t[5].split('/'))
			if ga[0] < ga[1]:
				data[(ga,gs)].append(l) #t[1])
				back.append((ga,gs))
			else:
				data[tuple(reversed(ga)),tuple(-s for s in reversed(gs))].append(l) #t[1])
				back.append((tuple(reversed(ga)),tuple(-s for s in reversed(gs))))

	f.close()
	return (data,back)

(data1,back1) = loadData(arguments["list1"])
(data2,back2) = loadData(arguments["list2"])

for s in back1:
	assert s in data1
	if s in data2:
#for s in set(data1).intersection(data2):
		g1 = data1[s].pop()
		for g2 in data2[s]:
		#for (g1,g2) in itertools.product(data1[s], data2[s]):
			print "%s\t%s" % (g1,g2)

