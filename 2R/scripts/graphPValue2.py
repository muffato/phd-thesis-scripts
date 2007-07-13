#! /users/ldog/muffato/python -OO

# INITIALISATION

# Librairies
import sys
import math

import utils.myMaths
import utils.myTools
import utils.walktrap

# Calcule la proba d'avoir l elements parmi ll, sachant qu'on en devrait avoir une proportion pi
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
	x = probaLog(p, para[c1][c2], nbPara/2)
	if para[c1][c2] > p*nbPara:
		x *= -1
	pvalues[c1][c2] = x
	#print c1, c2, x

#sys.exit(0)

print "graph {"
walktrap = utils.walktrap.WalktrapLauncher()

# Les liens entre chromosomes
#for (c1,c2) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.UpperMatrix):
for (c1,c2) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
	w = abs(pvalues[c1][c2])
	if w >= options["seuilPValue"]:
		#print "%d.%d -- %d.%d  [weight=%f,w=%f]" % (c1,c2, c2,c1, w,w)
		print "%d -- %d  [weight=%f,w=%f]" % (c1,c2, w,w)
	#walktrap.addEdge("%d.%d" % (c1,c2), "%d.%d" % (c2,c1), w)
	

# Les liens intra
for c1 in lstChr:
	continue
	for (c2,c3) in utils.myTools.myMatrixIterator(lstChr, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		
		#p = float(para[c1][c2] * para[c1][c3]) / float(nbPara*(nbPara-1))
		p = float(sum(para[c1].values())) * float(sum(para[c1].values())-1) / float(nbPara*(nbPara-1))

		if p > 0:
			w = abs(probaLog(p, min(para[c2][c3],min(para[c1][c2],para[c1][c3])), nbPara/2))
			#print >> sys.stderr, p, min(para[c2][c3],min(para[c1][c2],para[c1][c3])), nbPara/2, w
		else:
			w = 0

		#w = abs(pvalues[c2][c3])
		#if w >= options["seuilPValue"]:
		w = 1
		if w > 0:
			print "%d.%d -- %d.%d  [weight=%f,w=%f]" % (c1,c2, c1,c3, w,w)
			walktrap.addEdge("%d.%d" % (c1,c2), "%d.%d" % (c1,c3), w)
print "}"

sys.exit(0)
walktrap.doWalktrap()

for (nodes,cuts,_,dend) in walktrap.res:
	print >> sys.stderr, "Composante connexe: %d elements" % len(nodes)
	res = [(alpha,relevance,dend.cut(alpha)) for (alpha,relevance) in cuts]
	for (alpha,relevance,(clusters,lonely)) in res:
		print >> sys.stderr, "> alpha=%f relevance=%f clusters=%d size=%d lonely=%d sizes=%s" % \
			(alpha,relevance,len(clusters),sum([len(c) for c in clusters]),len(lonely),utils.myMaths.myStats([len(c) for c in clusters]))
		print >> sys.stderr, clusters

