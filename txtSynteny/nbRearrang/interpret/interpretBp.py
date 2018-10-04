#!/usr/bin/env python2

__doc__ = """
	Categorise les rearrangements en
		- insertion d'un gene
		- deplacement d'un gene unique avec conservation de l'orientation
"""

import sys
import itertools
import collections

import utils.myFile
import utils.myTools


arguments = utils.myTools.checkArgs( [("bpFile",file)], [], __doc__)

allGenes = set()

lstPresent = set()
lstAbsent = set()

def revGene(g):
	return (g[0],-g[1])

def revInt(i):
	return (revGene(i[1]), revGene(i[0]))

data = []
f = utils.myFile.openFile(arguments["bpFile"], "r")
for l in f:
	#l = l.replace("?", "+")
	t = l[:-1].split("\t")
	for (i,x) in enumerate(t):
		if "/" in x:
			break
	pattern = t[:i]
	#if "-" not in pattern:
	#	continue
	orient = [int(x) for x in t[i].split("/")]
	j = (len(t)-i-1)/2
	
	gene1 = (tuple(t[i+1:i+j+1]), orient[0])
	gene2 = (tuple(t[i+j+1:]), orient[1])
	allGenes.add(gene1)
	allGenes.add(gene2)

	intA = (gene1, gene2)
	intB = revInt(intA)

	data.append((pattern, intA))

	for (k,x) in enumerate(pattern):
		if x == "?":
			continue
		(lstPresent if x == "+" else lstAbsent).add( (intA, k) )
		(lstPresent if x == "+" else lstAbsent).add( (intB, k) )

f.close()

print >> sys.stderr, len(allGenes), len(data)

dicPresent = collections.defaultdict(dict)
for ((x,y),esp) in lstPresent:
	assert x not in dicPresent[esp], (esp,x,y,dicPresent[esp][x])
	dicPresent[esp][x] = y
print >> sys.stderr, len(lstPresent)/2

dicAbsent = collections.defaultdict(lambda : collections.defaultdict(set))
for ((x,y),esp) in lstAbsent:
	dicAbsent[esp][x].add(y)
print >> sys.stderr, len(lstAbsent)/2

lstEsp = sorted(set(dicPresent).union(dicAbsent))[1:]

type1 = set()
for a in dicPresent[0]:
	b = dicPresent[0][a]
	for esp in lstEsp:
		if a not in dicPresent[esp]:
			continue
		c = dicPresent[esp][a]
		if ((c,b),esp) in lstPresent:
			assert ((a,c),0) not in lstPresent
			assert ((c,b),0) not in lstPresent
			assert ((a,b),esp) not in lstPresent
			type1.add( (a,b,c,esp) )
			#print "[1]", a
			#print "[1]", b
			#print "[1]", c
			#print "[1]", esp
			#print 

print >> sys.stderr, "genes inseres seuls", len(type1)/2

type2 = set()
for esp in lstEsp:
	for a in dicPresent[esp]:
		b = dicPresent[esp][a]
		if a not in dicPresent[0]:
			continue
		c = dicPresent[0][a]
		if ((c,b),0) in lstPresent:
			assert ((a,c),esp) not in lstPresent
			assert ((c,b),esp) not in lstPresent
			assert ((a,b),0) not in lstPresent
			type2.add( (a,b,c,esp) )
			#print "[2]", a
			#print "[2]", b
			#print "[2]", c
			#print "[2]", esp
			#print 

print >> sys.stderr, "genes partis seuls", len(type2)/2

print >> sys.stderr, "intersection (genes deplaces seuls)", len(set((c,l) for (a,b,c,l) in type1).intersection((c,l) for (a,b,c,l) in type2))/2

changeToPresent = set()
changeToAbsent = set()
for (a,b,c,l) in type1:
	changeToPresent.add( ((a,b),esp) )
	changeToAbsent.add( ((a,c),esp) )
	changeToAbsent.add( ((c,b),esp) )

	#changeToAbsent.add( ((a,b),0) )
	#changeToPresent.add( ((a,c),0) )
	#changeToPresent.add( ((c,b),0) )

for (a,b,c,l) in type2:
	#changeToPresent.add( ((a,b),0) )
	#changeToAbsent.add( ((a,c),0) )
	#changeToAbsent.add( ((c,b),0) )

	changeToAbsent.add( ((a,b),esp) )
	changeToPresent.add( ((a,c),esp) )
	changeToPresent.add( ((c,b),esp) )

for (pattern,intA) in data:
	intB = revInt(intA)
	newPattern = ["+" if (intA,i) in changeToPresent else "-" if (intA,i) in changeToAbsent else x for (i,x) in enumerate(pattern)]
	print utils.myFile.myTSV.printLine(newPattern + ["%d/%d" % (intA[0][1],intA[1][1])] + list(intA[0][0]) + list(intA[1][0]))



