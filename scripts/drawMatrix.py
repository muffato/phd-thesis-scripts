#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine la matrice des genes orthologues entre deux genomes.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myDiags
import utils.myPsOutput


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("reverse",bool,False), ("scaleY",bool,False),\
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.loadGenome(options["orthologuesList"])
else:
	genesAnc = None
try:
	colors = utils.myGenomes.loadGenome(options["colorFile"])
except Exception:
	colors = None

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)


table12 = genome1.buildOrthosTable(chr1, genome2, chr2, options["includeGaps"], genesAnc)

# Matrice

print >> sys.stderr, "Affichage ",

utils.myPsOutput.printPsHeader()
if len(options["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, 21,29.7, options["backgroundColor"], options["backgroundColor"])
sys.stderr.write('.')

# Initialisations
table21 = genome2.buildOrthosTable(chr2, genome1, chr1, options["includeGaps"], genesAnc)
nb = sum([len(table12[c]) for c in table12])
scaleX = 19. / float(nb)
if options["scaleY"]:
	scaleY = 19. / float(sum([len(table21[c]) for c in table21]))
else:
	scaleY = scaleX
if options["pointSize"] < 0:
	dp = scaleX
else:
	dp = options["pointSize"]
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

lstNum1 = prepareGenome(table12, lambda x: utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, float(sum([len(table21[c]) for c in table21]))*scaleY, options["penColor"]))
sys.stderr.write('.')
lstNum2 = prepareGenome(table21, lambda y: utils.myPsOutput.drawLine(1, 1 + y*scaleY, 19, 0, options["penColor"]))
sys.stderr.write('.')

print "0 setlinewidth"

for c1 in table12:
	for (i1,t) in table12[c1]:
		xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
		for (c2,i2) in t:

			if colors == None:
				coul = options["defaultColor"]
			else:
				tmp = colors.getPosition(genome1.lstGenes[c1][i1])
				tmp.update( colors.getPosition(genome2.lstGenes[c2][i2]) )
				if len(tmp) == 0:
					coul = options["defaultColor"]
				else:
					coul = tmp.pop()[0]
			
			yy = 1 + lstNum2[(c2,i2)]*scaleY
			utils.myPsOutput.drawBox( xx, yy, dp, dp, coul, coul)

utils.myPsOutput.printPsFooter()
print >> sys.stderr, " OK"

