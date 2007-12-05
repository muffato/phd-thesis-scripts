#! /users/ldog/muffato/python -OO

import sys
import math
import operator
import utils.myMaths
import utils.myTools
import utils.myGenomes

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "paralogsList"], \
	[], \
	__doc__
)

genome = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])
paralogues = utils.myGenomes.Genome(noms_fichiers["paralogsList"])

para = utils.myTools.defaultdict(lambda : utils.myTools.defaultdict(int))
nbPara = 0
for g in paralogues:
	tg = genome.getPosition(g.names)
	for (c1,c2) in utils.myTools.myIterator.tupleOnStrictUpperList([c for (c,_) in tg if type(c) == int]):
		para[c1][c2] += 1
		para[c2][c1] += 1
		nbPara += 2
	#arcs.update( utils.myTools.myIterator.tupleOnStrictUpperList(list(tg)) )


# Arguments
#(noms_fichiers, options) = utils.myTools.checkArgs( [], [("seuilPValue",float,5)], "Lit une liste de paralogues (comme generee par convAncGenes.py et calcule une table de p-values")

#para = utils.myTools.defaultdict(lambda : utils.myTools.defaultdict(int))
#nbPara = 0
#for l in sys.stdin:
#	c = l.split()
#	c0 = int(c[0])
#	c2 = int(c[2])
#	para[c0][c2] += 1
#	para[c2][c0] += 1
#	nbPara += 2


lstChr = [x for x in para if len(para[x]) > 1]

pvalues = utils.myTools.defaultdict(dict)

for (c1,c2) in utils.myTools.myIterator.tupleOnWholeList(lstChr):
	p = float(sum(para[c1].values()))
	if c1 == c2:
		p *= float(sum(para[c1].values())-1)
	else:
		p *= float(sum(para[c2].values()))
	p /= float(nbPara*(nbPara-1))
	x = utils.myMaths.probaLog(p, para[c1][c2], nbPara/2)
	#if para[c1][c2] > p*nbPara:
	#	x *= -1
	pvalues[c1][c2] = abs(x)


res = {}
for n in xrange(3,8):
	for t in utils.myTools.myIterator.buildSubsets(lstChr, n):
		t = tuple(t)
		s = 0.
		nb = 0.
		for (x1,x2) in utils.myTools.myIterator.tupleOnStrictUpperList(t):
			s += pvalues[x1][x2]
			nb += 1
		res[t] = (s/nb, s)
	print >> sys.stderr, n, "ok"
print >> sys.stderr, len(res)

res2 = res.items()
res2.sort(reverse=True, key=operator.itemgetter(1))
for x in res2:
	print x

