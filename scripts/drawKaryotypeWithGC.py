#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import itertools
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPsOutput


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["studiedGenome", "referenceGenome", "GCGenesNames", "GCPercent"], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), \
	("GCorthologues",str,""), ("GCcolumn",int,0), ("GCaxisMin",float,25), ("GCaxisMax",float,85), ("GCsmoothing",int,10), \
	("reverse",bool,False), ("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("landscape",bool,False), ("showText",bool,True), ("drawBorder",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)


# Le genome avec les couleurs
genome2 = utils.myGenomes.Genome(noms_fichiers["referenceGenome"])

# Si on a utilise /dev/null, c'est que le caryotype est donne sous un autre format
if len(genome2.dicGenes) == 0:

	f = utils.myTools.myOpenFile(noms_fichiers["studiedGenome"], "r")
	table12 = {}
	for (i,l) in enumerate(f):
		nbChr = i+1
		chrom = []
		t = l.replace('\n', '').split(";")
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

# Sinon, procedure normale: on charge le genome avec les orthologues
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

(largeur,hauteur) = utils.myPsOutput.printPsHeader(landscape = True)

if len(options["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, largeur,hauteur, options["backgroundColor"], options["backgroundColor"])

# On dessine
if options["dx"] > 0:
	dx = options["dx"]
else:
	dx = (largeur-4.) / (5./2.*len(chr1) + 1.*(len(chr1)-1.))
if options["dy"] > 0:
	dy = options["dy"]
else:
	dy = (hauteur-4.) / float(max([len(x) for x in table12.values()]))
y0 = 1.

leniter = utils.myTools.myIterator.leniter
drawBox = utils.myPsOutput.drawBox


# Chargement du fichier avec les taux de GC
gcP = utils.myTools.myOpenFile(noms_fichiers["GCPercent"], "r")
gcN = utils.myTools.myOpenFile(noms_fichiers["GCGenesNames"], "r")
if len(options["GCorthologues"]) == 0:
	gcO = None
else:
	gcO = utils.myGenomes.Genome(options["GCorthologues"])
dicGC = {}
for (ligne1,ligne2) in itertools.izip(gcP,gcN):
	try:
		gc = float(ligne1.split()[options["GCcolumn"]])
	except ValueError:
		continue
	names = ligne2.split()
	if gcO != None:
		names = utils.myMaths.flatten( [gcO.lstGenes[c][i].names for (c,i) in gcO.getPosition(names)] )
	for (c,i) in genome1.getPosition(names):
		dicGC[(c,i)] = gc
gcP.close()
gcN.close()

dicGC2 = {}
for (c,i) in dicGC:
	tmp = [dicGC[(c,j)] for j in xrange(i-options["GCsmoothing"],i+1+options["GCsmoothing"]) if (c,j) in dicGC]
	dicGC2[(c,i)] = utils.myMaths.mean(tmp)


count = 0.
countNB = 0.
xx = 1
for c in chr1:
	def printBorder():
		print "newpath"
		print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,len(table12[c])*dy-dx)
		print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(table12[c])*dy-dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,-len(table12[c])*dy+dx)
		print "closepath"

	if options["roundedChr"]:
		print "initclip"
		printBorder()
		print "clip",

	def trans( (_,val) ):
		if len(val) == 0:
			return options["defaultColor"]
		else:
			return val[0][0]
	
	utils.myPsOutput.drawLine(xx+3./2.*dx, y0+0.75, dx, 0, options["penColor"])
	utils.myPsOutput.drawLine(xx+3./2.*dx, y0+0.7, 0, 0.1, options["penColor"])
	utils.myPsOutput.drawLine(xx+5./2.*dx, y0+0.7, 0, 0.1, options["penColor"])
	
	y = y0 + 1
	#lastGC = None
	for (col,items) in itertools.groupby(table12[c], key=trans):
		items = list(items)
		hauteur = len(items) * dy
		drawBox(xx, y, dx, hauteur, col, col)
		for (i,_) in items:
			GC = dicGC2.get( (c,i), None )
			if GC != None:
				GC = (GC-options["GCaxisMin"]) / (options["GCaxisMax"] - options["GCaxisMin"])
			if GC != None:
				drawBox(xx+3./2.*dx, y, dx*GC, dy, options["penColor"], options["penColor"])
			#if lastGC != None:
			#	utils.myPsOutput.drawLine(xx+3./2.*dx + dx*lastGC/100., y-dy, dx*(GC-lastGC)/100., dy, options["penColor"])
			
			#if (GC != None) and (lastGC != None):
			#	count += abs(lastGC-GC)
			#	countNB += 1.
			#if GC != None:
			#	lastGC = GC
			y += dy
	
	utils.myPsOutput.drawLine(xx+3./2.*dx, y+0.25, dx, 0, options["penColor"])
	utils.myPsOutput.drawLine(xx+3./2.*dx, y+0.2, 0, 0.1, options["penColor"])
	utils.myPsOutput.drawLine(xx+5./2.*dx, y+0.2, 0, 0.1, options["penColor"])
	
	print "initclip"
	if options["drawBorder"]:
		printBorder()
		utils.myPsOutput.setColor(options["penColor"], "color")
		print "stroke"

	if options["showText"]:
		utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
	
	xx += (7./2.*dx)

utils.myPsOutput.printPsFooter()
print >> sys.stderr, "OK" #, count/countNB


