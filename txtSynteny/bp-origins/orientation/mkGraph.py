#!/usr/bin/env python2

import sys
import math
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree",file)], [], "Cree un fichier xmgrace montrant la phylogenie")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])

resCoor = collections.defaultdict(list)
resConv = collections.defaultdict(list)
resDiv = collections.defaultdict(list)
for l in sys.stdin:
	#t = l[:-1].split("\t")
	t = l[:-1].split()

	#t2 = tuple([float(x) for x in t[1:]])
	#resCoor[t[0]] = (t2[4],0)
	#resConv[t[0]] = (t2[5],0)
	#resDiv[t[0]] = (t2[6],0)
	#resCoor[t[0]] = t2[:2]
	#resConv[t[0]] = t2[2:4]
	#resDiv[t[0]] = t2[4:6]

	if t[3] == t[4]:
		resCoor[t[0]].append(float(t[7]))
	elif t[3] == "-1":
		resDiv[t[0]].append(float(t[7]))
	else:
		resConv[t[0]].append(float(t[7]))
	
	#resCoor[t[1]].append(float(t[8]))
	#resConv[t[1]].append(float(t[9]))
	#resDiv[t[1]].append(float(t[10]))
	#resCoor[t[0]].append(float(t[6]))
	#resConv[t[0]].append(float(t[7]))
	#resDiv[t[0]].append(float(t[8]))

def do1(anc, lst):
	#return
	if anc not in lst:
		return
	#lst[anc] = [math.log(x/(bl[anc]/100.)+1e-3,2) for x in lst[anc]]
	#lst[anc] = (len([x for x in lst[anc] if x >= 2])/(bl[anc]/100.),0)
	#lst[anc] = [x/(bl[anc]/100.) for x in lst[anc]]
	m = utils.myMaths.myStats.mean(lst[anc])
	s = utils.myMaths.myStats.stddev(lst[anc])
	#print >> sys.stderr, anc, lst[anc], m, s, len(lst[anc]), math.sqrt(len(lst[anc])), (float(s)*1.96)/float(math.sqrt(len(lst[anc])))
	lst[anc] = (m, s*1.96/math.sqrt(len(lst[anc])))

n = 0
def do2(anc, lst, color):
	if anc not in lst:
		return
	if anc not in phylTree.parent:
		return
	par = phylTree.parent[anc][0]
	while (par not in lst) and (par in phylTree.parent):
		par = phylTree.parent[par][0]
	if par not in lst:
		return
	global n
	#print >> sys.stderr, anc, ">", par
	f = open('tmp/%d' % n, 'w')
	#print anc, par, lst[anc], lst[par]
	print >> f, phylTree.ages[anc], lst[anc][0], lst[anc][1]
	print >> f, phylTree.ages[par], lst[par][0], lst[par][1]
	#print >> f, phylTree.ages[anc], lst[anc][0], '"%s"' % anc
	#print >> f, phylTree.ages[par], lst[par][0], '"%s"' % par
	f.close()
	#print 'READ XY "tmp/%d"' % n
	print 'READ XYDY "tmp/%d"' % n
	print 's%d line color %d' % (n, color)
	print 's%d linewidth 2' % n
	print 's%d errorbar color %d' % (n, color)
	#print 's%d symbol 1' % n
	#print 's%d symbol size 0.3' % n
	#print 's%d symbol color %d' % (n, color)
	#print 's%d symbol fill color %d' % (n, color)
	#print 's%d symbol fill pattern 1' % n
	n += 1

bl = {}
def doBranchLength(node):
	s = 0.
	if node in phylTree.items:
		for (e,l) in phylTree.items[node]:
			s += l
			s += doBranchLength(e)
	bl[node] = s
	return s
doBranchLength(phylTree.root)
#print >> sys.stderr, bl

for anc in phylTree.listAncestr:
	do1(anc, resCoor)
	do1(anc, resConv)
	do1(anc, resDiv)

lstAnc = phylTree.listAncestr
lstAnc = [anc for anc in phylTree.listAncestr if (phylTree.dicParents[anc]["Euteleostomi"] == "Euteleostomi") and (anc != "Euteleostomi")]
#lstAnc = [anc for anc in lstAnc if len(phylTree.species[anc])>2]

for anc in lstAnc:
	do2(anc, resCoor, 1)
	do2(anc, resConv, 2)
	do2(anc, resDiv, 3)


