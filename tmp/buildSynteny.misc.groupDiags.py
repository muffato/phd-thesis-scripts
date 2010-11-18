#!/usr/bin/env python2

__doc__ = """
	Groupe toutes les paires de genes modernes orthologues et les paires de genes ancestrales correspondantes
"""


import sys
import itertools
import collections

import utils.myFile
import utils.myTools

arguments = utils.myTools.checkArgs( [("all-pairwise",file)], [], __doc__ )

#2       Erinaceus europaeus     GeneScaffold_7817       ENSEEUG00000006652 ENSEEUG00000006783   Xenopus tropicalis      scaffold_1080   ENSXETG00000014499 ENSXETG00000014497 1 1     Insectivora=4452 15629|Boreoeutheria=6591 22869|Theria=5731 19223|Eutheria=6392 21882|Mammalia=5627 18530|Amniota=5618 18305|Tetrapoda=5481 17639|Laurasiatheria=6061 21408

def trname(s):
	try:
		return int(s)
	except:
		return s

combin = utils.myTools.myCombinator()
f = utils.myFile.openFile(arguments["all-pairwise"], "r")
for l in f:
	t = l[:-1].split("\t")
	strands = [int(x) for x in t[7].split()]
	data = []
	data.append((t[1],t[3].split()))
	data.append((t[4],t[6].split()))
	for l in t[8].split("|"):
		(anc,_,l) = l.partition("=")
		data.append((anc,l.split()))
	names = [x[0] for x in data]
	lg = [[trname(y) for y in x[1]] for x in data]
	assert len(set(len(x) for x in lg)) == 1
	for ((s1,l1),(s2,l2)) in utils.myTools.myIterator.slidingTuple(zip(strands,zip(*lg))):
		tmp = []
		for (n,g1,g2) in itertools.izip(names, l1, l2):
			if g1 != g2:
				tmp.append( (n,g1,g2,s1,s2) )
				tmp.append( (n,g2,g1,-s2,-s1) )
		combin.addLink(tmp)
f.close()

for grp in combin:
	lst = [x for x in grp if x[1] < x[2]]
	if int in [type(x[1]) for x in lst]:
		print [x for x in grp if x[1] < x[2]]

