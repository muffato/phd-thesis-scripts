#!/usr/bin/env python2

import sys
import operator
import collections

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("boreoGenome",file), ("intervList",utils.myTools.FileList(1))], [("filterFirstSet",bool,True)], "")

boreoGenome = utils.myGenomes.Genome(arguments["boreoGenome"])

def loadIntervals(s):

	intervOK = set()
	intervNO = set()

	f = utils.myFile.openFile(s, "r")
	for line in f:
		interv = tuple(line.split()[1:3])
		if line.startswith("ANC-ADJ=NO") or line.startswith("ANC-ADJ=EXTR"):
		# Au cas ou un intervalle pourrait etre vu a la fois en OK et en NO
		# Dans ce cas, c'est le OK qui l'emporte
		#if line.startswith("anc-OK"):
		#	intervNO.discard(interv)
		#	intervOK.add(interv)
		#elif line.startswith("anc-NO") and (interv not in intervOK):
			intervNO.add(interv)
	f.close()
	return intervNO

comb = utils.myTools.myCombinator()
def addLinks(key, interv, filt=False):
	for (i,(xa,xb)) in enumerate(interv):
		(ca,ia) = boreoGenome.dicGenes[xa]
		(cb,ib) = boreoGenome.dicGenes[xb]
		assert ca == cb
		l = [(ca,x,boreoGenome.lstGenes[ca][x].names[0]) for x in (xrange(ia,ib+1) if ia < ib else xrange(ib, ia+1))]
		li = list(utils.myTools.myIterator.slidingTuple(l))
		added = False
		for x in li:
			if (not filt) or (x in count):
				count[x] += 1
				added = True
		if added:
			comb.addLink([key + str(i+1)] + li)

# Chargement de chaque jeu d'intervalles
count = collections.defaultdict(int)
for x in arguments["intervList"][1:]:
	interv = loadIntervals(x)
	addLinks(x, interv)
interv = loadIntervals(arguments["intervList"][0])
addLinks(arguments["intervList"][0], interv, filt=arguments["filterFirstSet"])
print >> sys.stderr, "score:", utils.myMaths.myStats.txtSummary(count.values())

print >> sys.stderr, len(list(comb)), "groups"
for grp in comb:

	# Intervalles tries selon leurs positions
	grp.sort()

	# On rajoute le score de chaque intervalle
	grp = [(count[x],x) for x in grp if isinstance(x, tuple)]
	assert len(grp) > 0
	assert min(grp)[0] > 0
	
	# Meilleur score
	bestS = max(grp)[0]

	# Les positions ou le maximum est atteint
	li = set(j for j in xrange(len(grp)) if grp[j][0] == bestS)
	
	# La region dans laquelle le rearrangement se trouve	
	print "REGION",
	for x in grp:
		print x[1][0][2], "[%d]" % x[0],
	print grp[-1][1][1][2]

	# On cherche les positions uniques
	for j in li:
		if ((j-1) not in li) and ((j+1) not in li):
			x = grp[j]
			print "BP", x[1][0][2], "[%d]" % x[0], x[1][1][2]

