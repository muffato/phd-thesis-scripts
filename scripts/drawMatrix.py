#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Dessine la matrice des genes orthologues entre deux genomes.
"""

import sys
import utils.myTools
import utils.myGenomes
import utils.myPsOutput

# Arguments
arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("referenceGenome",file)], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), ("includeNones",bool,False), \
	("reverse",bool,False), ("scaleY",bool,False),\
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"])
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])
if arguments["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if arguments["orthologuesList"] != "":
	genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"])
else:
	genesAnc = None
if arguments["colorFile"] != "":
	colors = utils.myGenomes.Genome(arguments["colorFile"])
else:
	colors = None

chr1 = []
chr2 = []
chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Chromosome])
chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Chromosome])
if arguments["includeScaffolds"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Scaffold])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Scaffold])
if arguments["includeRandoms"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Random])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Random])
if arguments["includeNones"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.None])
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.None])


table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)

# Matrice

print >> sys.stderr, "Affichage ",

utils.myPsOutput.printPsHeader()
if len(arguments["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, 21,29.7, arguments["backgroundColor"], arguments["backgroundColor"])
sys.stderr.write('.')

# Initialisations
table21 = genome2.buildOrthosTable(chr2, genome1, chr1, arguments["includeGaps"], genesAnc)
nb = sum([len(table12[c]) for c in table12])
scaleX = 19. / float(nb)
if arguments["scaleY"]:
	scaleY = 19. / float(sum([len(table21[c]) for c in table21]))
else:
	scaleY = scaleX
if arguments["pointSize"] < 0:
	dp = scaleX
else:
	dp = arguments["pointSize"]
sys.stderr.write('.')

def prepareGenome(dicOrthos, func):
	i = 0
	y = 0
	lstNum = {}
	func(y)
	for c in sorted(dicOrthos):
		y += len(dicOrthos[c])
		func(y)
		for (gene,_) in dicOrthos[c]:
			lstNum[(c,gene)] = i
			i += 1
	return lstNum

lstNum1 = prepareGenome(table12, lambda x: utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, float(sum([len(table21[c]) for c in table21]))*scaleY, arguments["penColor"]))
sys.stderr.write('.')
lstNum2 = prepareGenome(table21, lambda y: utils.myPsOutput.drawLine(1, 1 + y*scaleY, 19, 0, arguments["penColor"]))
sys.stderr.write('.')

print "0 setlinewidth"

for c1 in table12:
	for (i1,t) in table12[c1]:
		xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
		for (c2,i2) in t:

			if colors == None:
				coul = arguments["defaultColor"]
			else:
				tmp = set(colors.getPosition(genome1.lstGenes[c1][i1].names + genome2.lstGenes[c2][i2].names))
				if genesAnc != None:
					for (c,i) in genesAnc.getPosition(genome1.lstGenes[c1][i1].names + genome2.lstGenes[c2][i2].names):
						tmp.update(colors.getPosition(genesAnc.lstGenes[c][i].names))

				if len(tmp) == 0:
					coul = arguments["defaultColor"]
				else:
					coul = tmp.pop()[0]
			
			yy = 1 + lstNum2[(c2,i2)]*scaleY
			utils.myPsOutput.drawBox( xx, yy, dp, dp, coul, coul)

utils.myPsOutput.printPsFooter()
print >> sys.stderr, " OK"

