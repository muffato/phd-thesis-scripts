#! /users/ldog/muffato/python

__doc__ = """
	Affiche la matrice des p-values (d'existence/d'absence de paralogues) +/- en log
"""

import sys
import utils.myMaths
import utils.myTools
import utils.myGenomes

# Arguments
arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("paralogsList",file)], \
	[("greater",bool,True), ("log",bool,False), ("format",str,"%g")], \
	__doc__
)

genome = utils.myGenomes.Genome(arguments["studiedGenome"])
paralogues = utils.myGenomes.Genome(arguments["paralogsList"])

print >> sys.stderr, "Computing paralogs ...",
para = utils.myTools.defaultdict(lambda : utils.myTools.defaultdict(int))
nbPara = 0
for g in paralogues:
	tg = [c for (c,_) in genome.getPosition(g.names) if type(c) == int]
	for (c1,c2) in utils.myTools.myIterator.tupleOnStrictUpperList(tg):
		para[c1][c2] += 1
		para[c2][c1] += 1
		nbPara += 2
print >> sys.stderr, nbPara/2, "OK"

print >> sys.stderr, "Computing p-values ...",
lstChr = [x for x in para if len(para[x]) > 1]
pvalues = utils.myTools.defaultdict(dict)
for (c1,c2) in utils.myTools.myIterator.tupleOnWholeList(lstChr):
	p = float(sum(para[c1].values()))
	if c1 == c2:
		p *= float(sum(para[c1].values())-1)
	else:
		p *= float(sum(para[c2].values()))
	p /= float(nbPara*(nbPara-1))
	pvalues[c1][c2] = utils.myMaths.binomPvalue(p, para[c1][c2], nbPara/2, arguments["greater"])
print >> sys.stderr, "OK"

if arguments["log"]:
	import math
	transform = lambda x: arguments["format"] % abs(math.log10(x))
else:
	transform = lambda x: arguments["format"] % x

print utils.myFile.myTSV.printLine(["+"] + lstChr)
for c1 in lstChr:
	print utils.myFile.myTSV.printLine([c1] + [transform(pvalues[c1][c2]) for c2 in lstChr])

