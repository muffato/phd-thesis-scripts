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
	["lstSegments"], \
	[("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("showText",bool,False), ("drawBorder",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("centromereColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)

# Chargement des fichiers

f = open(noms_fichiers["lstSegments"], "r")
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

xx = 1
for c in chr1:
	def printBorder():
		print "newpath"
		print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,len(table12[c])*dy-dx)
		print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(table12[c])*dy-dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,-len(table12[c])*dy+dx)
		print "closepath"

	if options["dy"] < 0:
		dy = (hauteur-4.) / float(len(table12[c]))

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
		utils.myPsOutput.drawBox(xx, y, dx, hauteur, col, col)
		y += hauteur
	
	print "initclip"
	if options["drawBorder"]:
		printBorder()
		utils.myPsOutput.setColor(options["penColor"], "color")
		print "stroke"

	if options["showText"]:
		utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])

	y = y0 + 1
	for (col,items) in itertools.groupby(table12[c], key=trans):
		hauteur = leniter(items) * dy
		if col == "*CENTROMERE*":
			if len(options["backgroundColor"]) == 0:
				options["backgroundColor"] = "white"
			utils.myPsOutput.drawBox(xx-dx/3., y, dx*5./3., hauteur, options["backgroundColor"], options["backgroundColor"])
			print options["centromereColor"], "color"
			print "newpath"
			print "%.5f %.5f 2cm moveto" % (xx,y)
			print "%.5f %.5f 2cm rlineto" % (dx,hauteur)
			print "%.5f %.5f 2cm rlineto" % (-dx,0)
			print "%.5f %.5f 2cm rlineto" % (dx,-hauteur)
			print "closepath"
			print "fill"
		y += hauteur
	
	xx += (5.*dx)/3.

utils.myPsOutput.printPsFooter()
print >> sys.stderr, "OK"


