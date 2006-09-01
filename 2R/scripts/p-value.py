#!/usr/bin/python2.4

# INITIALISATION

# Librairies
import sys
import os
import math

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


def proba(pi, l, ll):
	p = pow(pi, l) * pow(1.-pi, ll-l)
	for i in range(l):
		p *= float(ll-i)/float(i+1)
	return p

def probaLog(pi, l, ll):
	p = 0.
	for i in range(l):
		p += math.log10(ll-i) - math.log10(i+1)
	#print >> sys.stderr, p, pi, l, ll
	p += l*math.log10(pi) + (ll-l)*math.log10(1-pi)
	return p



# Arguments
(noms_fichiers, options) = myTools.checkArgs( [], [], "")

para = {}
nbPara = 0
for l in sys.stdin:
	c = l.split()
	if c[0] not in para:
		para[c[0]] = {}
	if c[2] not in para:
		para[c[2]] = {}
	para[c[0]][c[2]] = para[c[0]].get(c[2], 0) + 1
	para[c[2]][c[0]] = para[c[2]].get(c[0], 0) + 1
	nbPara += 2


lstChr = para.keys()
lstChr.sort()

for c1 in lstChr:
	print "\t%s" % c1,
print

for c1 in lstChr:
	print c1,
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
	print


sys.exit(0)

genomeOutgroup = myOrthos.AncestralGenome(noms_fichiers[0], True)
genomePostDup = myOrthos.AncestralGenome(noms_fichiers[1], True)

# On genere toutes les tetrades
tetrades = myTools.buildSubsets(genomePostDup.lstChr, 4)
#tetrades.extend( myTools.buildSubsets(genomePostDup.lstChr, 3) )
#tetrades.extend( myTools.buildSubsets(genomePostDup.lstChr, 5) )
#tetrades.extend( myTools.buildSubsets(genomePostDup.lstChr, 6) )
print >> sys.stderr, len(tetrades), "tetrades"

# Un petit test, la somme des probabilites doit faire 1
s = 0.
for i in range(1001):
	s += proba(4./25., i, 1000)
print >> sys.stderr, "test des p-values:", s

# On fait la liste des orthologues et le tableau des comptes
count = {}
countChr = {}
totalCount = dict([ (k,0) for k in genomePostDup.lstChr ])
for c in genomeOutgroup.lstChr:
	count[c] = dict([ (k,0) for k in genomePostDup.lstChr ])
	countChr[c] = 0
	orthos = set([])
	for x in genomeOutgroup.lstGenes[c]:
		for g in x:
			if g in genomePostDup.dicGenes:
				orthos.add(genomePostDup.dicGenes[g])
	for (col,_) in orthos:
		count[c][col] += 1
		totalCount[col] += 1
		countChr[c] += 1
nbOrthos = sum(totalCount.values())
scoreTet = []
for tet in tetrades:
	nb = 0
	for k in tet:
		nb += totalCount[k]
	scoreTet.append( float(nb)/float(nbOrthos) )
print >> sys.stderr, nbOrthos, "orthologues"
print >> sys.stderr, totalCount

# On genere les probas
lstValue = {}
for c in genomeOutgroup.lstChr:
	cc = count[c]
	l = len(genomeOutgroup.lstGenes[c])
	lst = []
	for i in range(len(tetrades)):
		nb = 0
		for k in tetrades[i]:
			nb += cc[k]
		p = math.log10(proba(scoreTet[i], nb, l))
		lst.append(p)
		print int(100.*p)/100., c,
		for k in tetrades[i]:
			print k,
		print
	lstValue[c] = lst
	print >> sys.stderr, c, ":", countChr[c], "orthologues, log10(p-value) (min/max/moyenne):", min(lst), max(lst), myMaths.moyenne(lst)

sys.exit(0)
hist = {}
for i in range(len(tetrades)):
	lst = [lstValue[c][i] for c in genomeOutgroup.lstChr]
	for c in genomeOutgroup.lstChr:
	#lst2 = [math.log10(x) for x in lst]
	#if len([x for x in lst2 if x > -1.8]) != 0:
	#	continue
	
		x = lstValue[c][i]
		y = int(100.*math.log10(x))/100.
		print y, c,
		for k in tetrades[i]:
			print k,
		print
	continue
	#print lst2
	print 
	continue
	x = max(lst)
	y = int(100.*math.log10(x))/100.
	#y = myMaths.moyenne([int(100.*math.log10(x))/100. for x in lst])
	print y
	continue
	for c in genomeOutgroup.lstChr:
		x = lstValue[c][i]
		y = int(10.*math.log10(x))/10.
		if y in hist:
			hist[y] += 1
		else:
			hist[y] = 1

for y in hist:
	#print y, hist[y]
	print y, math.log10(hist[y])
