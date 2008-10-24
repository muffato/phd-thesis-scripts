#! /users/ldog/muffato/python

import sys
import math
import operator
import itertools
import collections

import utils.myMaths
import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myPsOutput
import utils.myKaryoDrawer

arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("refGenome",file), ("speciesName",str)], [("ancGenesFile1",str,""),("syntenyFile1",str,""),("ancGenesFile2",str,""),("syntenyFile2",str,"")], "")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
refGenome = utils.myGenomes.Genome(arguments["refGenome"])

refcolors = [(0,0,127), (0,192,192), (0,192,0), (255,255,0), (242,148,0), (224,0,0)]
inter = utils.myMaths.myInterpolator.getMultDim(utils.myMaths.myInterpolator.oneDimCubic, range(len(refcolors)), refcolors)

dataSynt1 = {}
dataSynt2 = {}
for (c,l) in refGenome.lstGenes.iteritems():
	dataSynt1[c] = [None] * len(l)
	dataSynt2[c] = [None] * len(l)

def addAncestr(anc, i):
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile1"] % phylTree.fileName[anc])
	f = utils.myFile.myTSV.reader(arguments["syntenyFile1"] % phylTree.fileName[anc])
	for l in f.csvobject:
		lg = [int(x) for x in l[2].split()]
		for g in lg:
			for (c,p) in refGenome.getPosition(ancGenes.lstGenes[None][g].names):
				dataSynt1[c][p] = i
	f.file.close()
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile2"] % phylTree.fileName[anc])
	f = utils.myFile.myTSV.fileReader(arguments["syntenyFile2"] % phylTree.fileName[anc])
	for l in f.csvobject:
		lg = [int(x) for x in l[2].split()]
		for g in lg:
			for (c,p) in refGenome.getPosition(ancGenes.lstGenes[None][g].names):
				dataSynt2[c][p] = i
	f.file.close()


num = 0
name = arguments["speciesName"]
dicName = {num: name}
while name in phylTree.parent:
	num += 1
	dicName[num] = name
	(name,_) = phylTree.parent[name]
	addAncestr(name, num)
col = utils.myPsOutput.getCubicGradient(refcolors, num+1)

all = collections.defaultdict(int)

newdata = {}
for c in refGenome.lstChr:

	size = len(refGenome.lstGenes[c])
	dummy = []
	for i in xrange(size):
		print dataSynt1[c][i], dataSynt2[c][i]
		#dummy.append(dataSynt1[c][i]-dataSynt2[c][i])
		#if (dataSynt1[c][i]*dataSynt2[c][i]) == 0:
		#	if dataSynt1[c][i] != dataSynt2[c][i]:
		#		print >> sys.stderr, dataSynt1[c][i], dataSynt2[c][i]
		if dataSynt1[c][i] == dataSynt2[c][i]:
			dataSynt1[c][i] = dataSynt2[c][i] = None

	#tmp = [utils.myMaths.myStats.mean(dataSynt1[c][max(i-4,0):i+5]) for i in xrange(size)]
	tmp = dataSynt1[c]
	dataSynt1[c] = [None if x is None else inter(x*(len(refcolors)-1.)/num) for x in tmp]

	#tmp = [utils.myMaths.myStats.mean(dataSynt2[c][max(i-4,0):i+5]) for i in xrange(size)]
	tmp = dataSynt2[c]
	dataSynt2[c] = [None if x is None else inter(x*(len(refcolors)-1.)/num) for x in tmp]

	newdata[c] = (dataSynt1[c],dataSynt2[c])
	
	#newdata[c] = [col[i] for i in data[c]]
	#for (x,l) in itertools.groupby(enumerate(data[c]), lambda x: x[1]):
	#	if x >= 10:
	#		print >> sys.stderr, c, dicName[x], [refGenome.lstGenes[c][i].names[0] for (i,_) in l]
	for (x,y) in utils.myMaths.myStats.count(dummy).iteritems():
		all[x] += y
	#print >> sys.stderr, c, utils.myMaths.myStats.count(dummy)
newdata[None] = utils.myMaths.flatten([[c]*50 for c in col])

#for (x,y) in all.iteritems():
#	print x, y

print >> sys.stderr, "Affichage ...",
#utils.myKaryoDrawer.drawKaryo(newdata, arguments)
print >> sys.stderr, "OK"

