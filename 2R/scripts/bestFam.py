#! /users/ldog/muffato/python -OO

# INITIALISATION

# Librairies
import sys
import math
import operator
import utils.myTools

# Calcule la proba
def proba(pi, l, ll):
	p = pow(pi, l) * pow(1.-pi, ll-l)
	for i in range(l):
		p *= float(ll-i)/float(i+1)
	return p


# Calcule le log de la proba
def probaLog(pi, l, ll):
	p = l*math.log10(pi) + (ll-l)*math.log10(1-pi)
	for i in range(l):
		p += math.log10(ll-i) - math.log10(i+1)
	return p



# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( [], [("seuilPValue",float,5)], "Lit une liste de paralogues (comme generee par convAncGenes.py et calcule une table de p-values")

para = utils.myTools.defaultdict(lambda : utils.myTools.defaultdict(int))
nbPara = 0
for l in sys.stdin:
	c = l.split()
	c0 = int(c[0])
	c2 = int(c[2])
	para[c0][c2] += 1
	para[c2][c0] += 1
	nbPara += 2


lstChr = para.keys()

nm = max(lstChr)+1
pvalues = [[0.] * nm for i in xrange(nm)]

for (c1,c2) in utils.myTools.myIterator.tupleOnWholeList(lstChr):
	p = float(sum(para[c1].values()))
	if c1 == c2:
		p *= float(sum(para[c1].values())-1)
	else:
		p *= float(sum(para[c2].values()))
	p /= float(nbPara*(nbPara-1))
	x = probaLog(p, para[c1][c2], nbPara/2)
	#if para[c1][c2] > p*nbPara:
	#	x *= -1
	pvalues[c1][c2] = abs(x)


res = {}
for l in xrange(2,7):
	lst = [lstChr] * l
	for t in utils.myTools.myIterator.tupleOnManyLists(*lst):
		tt = frozenset(t)
		if len(tt) != l:
			continue
		if tt in res:
			continue
		s = 0.
		nb = 0.
		for (x1,x2) in utils.myTools.myIterator.tupleOnStrictUpperList(t):
			s += pvalues[x1][x2]
			nb += 1
		res[tt] = (s/nb, s, nb, t)
	print >> sys.stderr, l, "ok"

res2 = res.items()
res2.sort(reverse=True, key=operator.itemgetter(1))
for x in res:
	print x

