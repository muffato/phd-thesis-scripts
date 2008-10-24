#! /users/ldog/muffato/python

import sys
import math
import operator
import itertools

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myPsOutput
import utils.myKaryoDrawer

arguments = utils.myTools.checkArgs([("phylTree.conf",file), ("refGenome",file), ("speciesName",str)], [("ancGenesFile",str,""),("syntenyFile",str,"")], "")

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
refGenome = utils.myGenomes.Genome(arguments["refGenome"])

refcolors = [(0,0,127), (0,192,192), (0,192,0), (255,255,0), (242,148,0), (224,0,0)]
inter = utils.myMaths.myInterpolator.getMultDim(utils.myMaths.myInterpolator.oneDimCubic, range(len(refcolors)), refcolors)

dataAge = {}
dataSynt = {}
for (c,l) in refGenome.lstGenes.iteritems():
	dataAge[c] = [0] * len(l)
	dataSynt[c] = [0] * len(l)

def addAncestr(anc, i):
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	for g in ancGenes:
		for (c,p) in refGenome.getPosition(g.names):
			dataAge[c][p] = i
	f = utils.myFile.myTSV.fileReader(arguments["syntenyFile"] % phylTree.fileName[anc])
	for l in f.csvobject:
		lg = [int(x) for x in l[2].split()]
		for g in lg:
			for (c,p) in refGenome.getPosition(ancGenes.lstGenes[None][g].names):
				dataSynt[c][p] = i
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

#addAncestrGenes("Homo/Pan/Gorilla group", 1)
#addAncestrGenes("Boreoeutheria", 2)
#addAncestrSynteny("Boreoeutheria", 3)
#col = [0, (0,0,127), (0,192,192), (0,192,0)]


newdata = {}
for c in refGenome.lstChr:

	size = len(refGenome.lstGenes[c])

	diff = []
	for i in xrange(size):
		assert dataSynt[c][i] <= dataAge[c][i]
		#diff.append(dataAge[c][i]-dataSynt[c][i])
		diff.append(phylTree.ages[dicName[dataAge[c][i]]] - phylTree.ages[dicName[dataSynt[c][i]]])
		diff[-1] = (math.log(diff[-1]+1) * num) / math.log(phylTree.ages[phylTree.root]+1)

	tmp = [utils.myMaths.myStats.mean(diff[max(i-4,0):i+5]) for i in xrange(size)]
	#tmp = diff
	diff = [inter(x*(len(refcolors)-1.)/num) for x in tmp]

	tmp = [utils.myMaths.myStats.mean(dataAge[c][max(i-4,0):i+5]) for i in xrange(size)]
	#tmp = dataAge[c]
	dataAge[c] = [inter(x*(len(refcolors)-1.)/num) for x in tmp]

	tmp = [utils.myMaths.myStats.mean(dataSynt[c][max(i-4,0):i+5]) for i in xrange(size)]
	#tmp = dataSynt[c]
	dataSynt[c] = [inter(x*(len(refcolors)-1.)/num) for x in tmp]

	newdata[c] = diff
	#newdata[c] = (dataAge[c],dataSynt[c])
	
	#newdata[c] = [col[i] for i in data[c]]
	#for (x,l) in itertools.groupby(enumerate(data[c]), lambda x: x[1]):
	#	if x >= 10:
	#		print >> sys.stderr, c, dicName[x], [refGenome.lstGenes[c][i].names[0] for (i,_) in l]
	#print >> sys.stderr, c, utils.myMaths.myStats.count(data[c])
newdata[None] = utils.myMaths.flatten([[c]*50 for c in col])

print >> sys.stderr, "Affichage ...",
utils.myKaryoDrawer.drawKaryo(newdata, arguments)
print >> sys.stderr, "OK"

