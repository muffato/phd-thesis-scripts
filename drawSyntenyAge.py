#! /users/ldog/muffato/python

import sys
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

data = {}
for (c,l) in refGenome.lstGenes.iteritems():
	data[c] = [0] * len(l)

def addAncestrGenes(anc, i):
	for g in utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc]):
		for (c,p) in refGenome.getPosition(g.names):
			data[c][p] = i


def addAncestrSynteny(anc, i):
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])
	f = utils.myTools.tsvReader(arguments["syntenyFile"] % phylTree.fileName[anc])
	for l in f.csvobject:
		lg = [int(x) for x in l[2].split()]
		for g in lg:
			for (c,p) in refGenome.getPosition(ancGenes.lstGenes[None][g].names):
				data[c][p] = i
	f.file.close()


num = 0
name = arguments["speciesName"]
while name in phylTree.parent:
	num += 1
	(name,_) = phylTree.parent[name]
	addAncestrGenes(name, num)
col = utils.myPsOutput.getCubicGradient(refcolors, num+1)

#addAncestrGenes("Homo/Pan/Gorilla group", 1)
#addAncestrGenes("Boreoeutheria", 2)
#addAncestrSynteny("Boreoeutheria", 3)
#col = [0, (0,0,127), (0,192,192), (0,192,0)]

inter = utils.myMaths.myInterpolator.getMultDim(utils.myMaths.myInterpolator.oneDimCubic, range(len(refcolors)), refcolors)

newdata = {}
for c in refGenome.lstChr:
	tmp = [utils.myMaths.myStats.mean(data[c][i-4:i+5])*5./num for i in xrange(len(data[c]))]
	newdata[c] = [inter(x) for x in tmp]
	#newdata[c] = [col[i] for i in data[c]]
	#for (i,x) in enumerate(data[c]):
	#	if x == 0:
	#		print >> sys.stderr, refGenome.lstGenes[c][i]
	#print >> sys.stderr, c, utils.myMaths.myStats.count(data[c])

print >> sys.stderr, "Affichage ...",
utils.myKaryoDrawer.drawKaryo(newdata, arguments)
print >> sys.stderr, "OK"

