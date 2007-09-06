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
	para[c[0]][c[2]] += 1
	para[c[2]][c[0]] += 1
	nbPara += 2


lstChr = para.keys()
lstChr.sort()

for c1 in lstChr:
	print "\t%s" % c1,
print

pvalues = {}
for c1 in lstChr:
	print c1,
	pvalues[c1] = {}
	for c2 in lstChr:
		p = float(sum(para[c1].values()))
		if c1 == c2:
			p *= float(sum(para[c1].values())-1)
		else:
			p *= float(sum(para[c2].values()))
		p /= float(nbPara*(nbPara-1))
		x = probaLog(p, para[c1].get(c2, 0), nbPara/2)
		if para[c1].get(c2, 0) > p*nbPara:
			x *= -1
		print "\t%g" % x,
		pvalues[c1][c2] = x
		if abs(x) >= options["seuilPValue"] and c1 > c2:
			print "%s %s %f" % (c1,c2,abs(x))
	print


sys.exit(0)
s0 = [set(x) for x in lstChr]
ss = [s0]

while True:
	newS = []
	for x in lstChr:
		for s in ss[-1]:
			if x in s:
				continue
			#m = 0.
			for c in s:
				#m += pvalues[x][c]
				if pvalues[x][c] < options["seuilPValue"]:
					break
			#m /= len(s)
			#if m >= options["seuilPValue"]:
			else:
				s2 = s.union([x])
				if s2 not in newS:
					newS.append(s2)
	if len(newS) == 0:
		break
	ss.append(newS)

newSS = []
for i in range(len(ss)-1):
	newSS.append([x for x in ss[i] if True not in [y.issuperset(x) for y in ss[i+1]]])
newSS.append(ss[-1])

for ts in newSS:
	print ts
		

#print >> sys.stderr, pvalues
