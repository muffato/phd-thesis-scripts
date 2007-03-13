#! /users/ldog/muffato/python -OO

__doc__ = """
Dessine la matrice des genes orthologues entre deux genomes.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPsOutput


#############
# FONCTIONS #
#############

# Fabrique la liste des orthologues utilises pour le dessin
def buildOrthosTable(genome1, chr1, genome2, chr2):

	# Tous les orthologues entre pour les chromosomes OK du genome 1
	res = {}
	for c1 in chr1:
		res[c1] = {}
		for g1 in genome1.lstGenes[c1]:
			tmp = set([])
			for s in g1.names:
				if s in genome2.dicGenes:
					tmp.add(genome2.dicGenes[s])
				if options["orthologuesList"] == "":
					continue
				if s in genesAnc.dicGenes:
					(c,i) = genesAnc.dicGenes[s]
					for ss in genesAnc.lstGenes[c][i].names:
						if ss in genome2.dicGenes:
							tmp.add(genome2.dicGenes[ss])
			# On ne garde que les chromsomes OK du genome 2
			tmp = [(c,i) for (c,i) in tmp if c in chr2]
			# +/- includeGaps
			if (not options["includeGaps"]) and (len(tmp) == 0):
				continue
			res[c1][genome1.dicGenes[s][1]] = [(c,i) for (c,i) in tmp if c in chr2]
			

	return res


########
# MAIN #
########

# Arguments
modes = ["Matrix", "Karyotype", "OrthosChr"]
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps", bool, False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("output",str,modes), ("scaleY",bool,False), \
	("pointSize",float,-1), ("colorFile",str,""), ("defaultColor",str,"black"), ("penColor",str,"black"), ("minHomology",int,90)], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.loadGenome(noms_fichiers["studiedGenome"])
genome2 = utils.myGenomes.loadGenome(noms_fichiers["referenceGenome"])
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.AncestralGenome(options["orthologuesList"])
else:
	genesAnc = genome2
try:
	colors = utils.myGenomes.loadGenome(options["colorFile"])
except Exception:
	colors = options["defaultColor"]
mode = modes.index(options["output"])

# Les chromosomes a etudier
chr1 = genome1.lstChr
chr2 = genome2.lstChr
if options["includeScaffolds"]:
	chr1.extend(genome1.lstScaff)
	chr2.extend(genome2.lstScaff)
if options["includeRandoms"]:
	chr1.extend(genome1.lstRand)
	chr2.extend(genome2.lstRand)


table12 = buildOrthosTable(genome1, chr1, genome2, chr2)

# Matrice
if (mode == 0):

	print >> sys.stderr, "Affichage ",
	
	utils.myPsOutput.initColor()
	utils.myPsOutput.printPsHeader()
	sys.stderr.write('.')

	# Initialisations
	table21 = buildOrthosTable(genome2, chr2, genome1, chr1)
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
			for gene in dicOrthos[c]:
				lstNum[(c,gene)] = i
				i += 1
		return lstNum

	lstNum1 = prepareGenome(table12, lambda x: utils.myPsOutput.drawLine(1 + x*scaleX, 1, 0, float(sum([len(table21[c]) for c in table21]))*scaleY, options["penColor"]))
	sys.stderr.write('.')
	lstNum2 = prepareGenome(table21, lambda y: utils.myPsOutput.drawLine(1, 1 + y*scaleY, 19, 0, options["penColor"]))
	sys.stderr.write('.')

	for c1 in table12:
		for i1 in table12[c1]:
			
			xx = 1 + float(lstNum1[(c1,i1)]) * scaleX
			
			for (c2,i2) in table12[c1][i1]:

				if type(colors) == str:
					coul = colors
				else:
					for gt in genome1.lstGenes[c1][i1].names + genome2.lstGenes[c2][i2].names:
						if gt in colors.dicGenes:
							coul = colors.dicGenes[gt][0]
							break
					else:
						continue
				
				yy = 1 + float(lstNum2[(c2,i2)]) * scaleY
				coul = utils.myPsOutput.getColor(str(coul), options["defaultColor"])
				utils.myPsOutput.drawBox( xx, yy, dp, dp, coul, coul)

	utils.myPsOutput.printPsFooter()
	print >> sys.stderr, " OK"


# Caryotype
elif (mode == 1):

	print >> sys.stderr, "Affichage ...",
	
	utils.myPsOutput.initColor()
	utils.myPsOutput.printPsHeader()

	# On dessine
	dx = (19.*3.)/(5.*len(chr1)-2.)
	dy = float(max([len(x) for x in table12.values()])) / 26.
	xx = 1
	y0 = 1.

	for c in chr1:
		
		if options["scaleY"]:
			dy = float(len(table12[c])) / 26.

		utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
		y = y0 + 1
		
		last = ""
		nb = 0
		for i in table12[c]:
			tmp = table12[c][i]
			if len(tmp) == 0:
				col = options["defaultColor"]
			else:
				col = utils.myPsOutput.getColor(str(tmp[0][0]), options["defaultColor"])
			if col == last:
				nb += 1
			else:
				if nb != 0:
					utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
					y += nb/dy
				last = col
				nb = 1
		if nb != 0:
			utils.myPsOutput.drawBox(xx, y, dx, nb/dy, last, last)
		xx += (5.*dx)/3.

	utils.myPsOutput.printPsFooter()
	print >> sys.stderr, "OK"


# Liste de chromosomes orthologues
elif (mode == 2):

	for c1 in chr1:

		s = str(c1)
		
		# On ne garde que les chromosomes orthologues majoritaires
		lst = [x for (x,_) in utils.myMaths.flatten(table12[c1].itervalues())]
		count = [(lst.count(x),x) for x in set(lst)]
		count.sort()
		count.reverse()
		nb = (len(lst)*options["minHomology"])/100
		for (n,c2) in count:
			s += "\t%s (%d)" % (c2,n)
			nb -= n
			if nb <= 0:
				break
		print s
