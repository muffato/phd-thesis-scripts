#! /users/ldog/muffato/python -OO

# INITIALISATION

# Librairies
import sys
import math

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


lstChr = sorted(para.keys())
pvalues = utils.myTools.defaultdict(dict)

for (c1,c2) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.WholeMatrix):
	p = float(sum(para[c1].values()))
	if c1 == c2:
		p *= float(sum(para[c1].values())-1)
	else:
		p *= float(sum(para[c2].values()))
	p /= float(nbPara*(nbPara-1))
	x = probaLog(p, para[c1].get(c2, 0), nbPara/2)
	if para[c1].get(c2, 0) > p*nbPara:
		x *= -1
	pvalues[c1][c2] = x
	#if abs(x) >= options["seuilPValue"] and c1 > c2:
	#if c1 > c2:
	#	print >> sys.stderr, "%s -- %s [weight=%f,w=%f]" % (c1,c2,abs(x),abs(x))
	#	print >> sys.stderr, "%s %s %f" % (c1,c2,abs(x))


print "graph {"

# Les liens entre chromosomes
for (c1,c2) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.UpperMatrix):
	w = abs(pvalues[c1][c2])
	if w < options["seuilPValue"]:
		continue
	print "%d.%d -- %d.%d  [weight=%f,w=%f]" % (c1,c2, c2,c1, w,w)

# Les liens intra
for c1 in lstChr:
	for (c2,c3) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		w = abs(pvalues[c2][c3])
		if w < options["seuilPValue"]:
			continue
		print "%d.%d -- %d.%d  [weight=%f,w=%f]" % (c1,c2, c1,c3, w,w)

print "}"
