#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("boreoGenome",file), ("intervList",utils.myTools.FileList(1))], [], "")

boreoGenome = utils.myGenomes.Genome(arguments["boreoGenome"])

def loadIntervals(s):
	r = []
	f = utils.myFile.openFile(s, "r")
	for line in f:
		r.append(tuple(line.split()))
	f.close()
	return r

interv = [loadIntervals(x) for x in arguments["intervList"]]

comb = utils.myTools.myCombinator()

def addLinks(src, interv):
	for (i,(xa,xb)) in enumerate(interv):
		(ca,ia) = boreoGenome.dicGenes[xa]
		(cb,ib) = boreoGenome.dicGenes[xb]
		assert ca == cb
		l = [g.names[0] for g in (boreoGenome.lstGenes[ca][ia:ib+1] if ia < ib else boreoGenome.lstGenes[ca][ib:ia+1])]
		l = list(utils.myTools.myIterator.slidingTuple(l))
		comb.addLink(l + [src + str(i+1)])

for (i,l) in enumerate(interv):
	addLinks("src%d-" % (i+1), l)

for grp in comb:
	lint = sorted([x for x in grp if isinstance(x, str)]) + ["-".join(x) for x in grp if isinstance(x, tuple)]
	print " ".join(lint)

