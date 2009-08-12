#! /users/ldog/muffato/python

import sys
import math
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree",file), ("orientFile",file)], [], "Cree un fichier xmgrace montrant la phylogenie")

dicEsp = {}
f = utils.myFile.openFile(arguments["orientFile"], "r")
for l in f:
	t = l[:-1].split("\t")
	#dicEsp[t[0]] = (float(t[1]), float(t[2]), float(t[3]))
	dicEsp[t[1]] = (int(t[3]), int(t[4]), int(t[5]))
f.close()

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree"])

resCoor = collections.defaultdict(list)
resConv = collections.defaultdict(list)
resDiv = collections.defaultdict(list)
for l in sys.stdin:
	t = l[:-1].split("\t")
	#resCoor[t[1]].append(float(t[8]))
	#resConv[t[1]].append(float(t[9]))
	#resDiv[t[1]].append(float(t[10]))
	#expCoor = min(dicEsp[t[1]][0], dicEsp[t[0]][0])
	#expConv = min(dicEsp[t[1]][1], dicEsp[t[0]][1])
	#expDiv = min(dicEsp[t[1]][2], dicEsp[t[0]][2])
	expCoor = min(dicEsp[t[2]][0], dicEsp[t[3]][0])
	expConv = min(dicEsp[t[2]][1], dicEsp[t[3]][1])
	expDiv = min(dicEsp[t[2]][2], dicEsp[t[3]][2])
	if expCoor*expConv*expDiv == 0:
		continue
	#print >> sys.stderr, expCoor, expConv, expDiv
	s = float(expCoor + expConv + expDiv)
	expCoor *= 100/s
	expConv *= 100/s
	expDiv *= 100/s
	#expConv = dicEsp[t[1]][1]
	#expDiv = dicEsp[t[1]][2]
	#resCoor[t[0]].append(float(t[6]))
	#resConv[t[0]].append(float(t[7]))
	#resDiv[t[0]].append(float(t[8]))
	#resCoor[t[0]].append(math.log(float(t[6])/expCoor + 1e-7, 2))
	#resConv[t[0]].append(math.log(float(t[7])/expConv + 1e-7, 2))
	#resDiv[t[0]].append(math.log(float(t[8])/expDiv + 1e-7, 2))
	resCoor[t[1]].append(math.log(float(t[8])/expCoor + 1e-7, 2))
	resConv[t[1]].append(math.log(float(t[9])/expConv + 1e-7, 2))
	resDiv[t[1]].append(math.log(float(t[10])/expDiv + 1e-7, 2))

def do1(anc, lst):
	if anc not in lst:
		return
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
	#par = phylTree.parent[anc][0]
	f = open('tmp/%d' % n, 'w')
	#print anc, par, lst[anc], lst[par]
	print >> f, phylTree.ages[anc], lst[anc][0], lst[anc][1]
	print >> f, phylTree.ages[par], lst[par][0], lst[par][1]
	f.close()
	print 'READ XYDY "tmp/%d"' % n
	print 's%d line color %d' % (n, color)
	print 's%d linewidth 2' % n
	print 's%d errorbar color %d' % (n, color)
	n += 1


for anc in phylTree.listAncestr:
	do1(anc, resCoor)
	do1(anc, resConv)
	do1(anc, resDiv)

for anc in phylTree.listAncestr:
	if phylTree.dicParents[anc]["Euteleostomi"] != "Euteleostomi":
		continue
	if anc == "Euteleostomi":
		continue
	do2(anc, resCoor, 1)
	do2(anc, resConv, 2)
	do2(anc, resDiv, 3)


