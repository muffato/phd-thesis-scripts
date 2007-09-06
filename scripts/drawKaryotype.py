#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import operator
import itertools
import utils.myGenomes
import utils.myTools
import utils.myPsOutput


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("reverse",bool,False), ("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("landscape",bool,False), ("showText",bool,True), ("drawBorder",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Le premier genome
genome2 = utils.myGenomes.Genome(noms_fichiers["referenceGenome"])

# Si on a utilise /dev/null, c'est que le caryotype est donne sous un autre format
if len(genome2.dicGenes) == 0:

	f = open(noms_fichiers["studiedGenome"], "r")
	table12 = {}
	for (i,l) in enumerate(f):
		nbChr = i+1
		chrom = []
		t = l[:-1].split(";")
		for s in t:
			if len(s) == 0:
				continue
			elif len(s) == 1:
				nbChr = s[0]
			c = s.split()
			chrom.extend( [[(c[0],None)]] *  int(10*float(c[1])) )
		table12[nbChr] = list(enumerate(chrom))
		table12[nbChr].reverse()

	chr1 = sorted(table12.keys())

else:
	genome1 = utils.myGenomes.Genome(noms_fichiers["studiedGenome"])


if options["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if options["orthologuesList"] != "":
	genesAnc = utils.myGenomes.Genome(options["orthologuesList"])
else:
	genesAnc = None

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

print >> sys.stderr, "Affichage ...",

(largeur,hauteur) = utils.myPsOutput.printPsHeader()

if len(options["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, largeur,hauteur, options["backgroundColor"], options["backgroundColor"])

# On dessine
if options["dx"] > 0:
	dx = options["dx"]
else:
	dx = ((largeur-2.)*3.)/(5.*len(chr1)-2.)
if options["dy"] > 0:
	dy = options["dy"]
else:
	dy = (hauteur-4.) / float(max([len(x) for x in table12.values()]))
y0 = 1.

leniter = utils.myTools.leniter
drawBox = utils.myPsOutput.drawBox

xx = 1
for c in chr1:
	def printBorder():
		print "newpath"
		print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,len(table12[c])*dy-dx)
		print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(table12[c])*dy-dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,-len(table12[c])*dy+dx)
		print "closepath"

	#if options["dy"] < 0:
	#	dy = (hauteur-4.) / float(len(table12[c]))

	if options["roundedChr"]:
		print "initclip"
		printBorder()
		print "clip"

	def trans( (_,val) ):
		if len(val) == 0:
			return options["defaultColor"]
		else:
			return val[0][0]
	
	y = y0 + 1
	for (col,items) in itertools.groupby(table12[c], key=trans):
		hauteur = leniter(items) * dy
		drawBox(xx, y, dx, hauteur, col, col)
		y += hauteur
	
	print "initclip"
	if options["drawBorder"]:
		printBorder()
		utils.myPsOutput.setColor(options["penColor"], "color")
		print "stroke"

	if options["showText"]:
		utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
	
	xx += (5.*dx)/3.

utils.myPsOutput.printPsFooter()
print >> sys.stderr, "OK"


