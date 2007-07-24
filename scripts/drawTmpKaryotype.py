#! /users/ldog/muffato/python -OO

__doc__ = """
	Dessine le karyotype d'un genome face a l'autre
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myPsOutput


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["lstSegments"], \
	[("dx",float,-1), ("scaleY",bool,False), ("roundedChr",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("centromereColor",str,"black"), ("backgroundColor",str,"")], \
	__doc__
)

# Chargement des fichiers

f = open(noms_fichiers["lstSegments"], "r")
nbChr = 0
table12 = {}
for l in f:
	chrom = []
	t = l[:-1].split(";")
	for s in t:
		if len(s) == 0:
			continue
		c = s.split()
		chrom.extend( [[(c[0],None)]] *  int(10*float(c[1])) )
	nbChr += 1
	table12[nbChr] = list(enumerate(chrom))

chr1 = sorted(table12.keys())

print >> sys.stderr, "Affichage ...",

utils.myPsOutput.printPsHeader()

if len(options["backgroundColor"]) > 0:
	utils.myPsOutput.drawBox(0,0, 21,29.7, options["backgroundColor"], options["backgroundColor"])

# On dessine
if options["dx"] > 0:
	dx = options["dx"]
else:
	dx = (19.*3.)/(5.*len(chr1)-2.)
dy = float(max([len(x) for x in table12.values()])) / 26.
y0 = 1.

if options["roundedChr"]:
	xx = 1
	print "newpath"
	for c in chr1:
		y = y0+1
		print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,len(table12[c])/dy-dx)
		print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(table12[c])/dy-dx/2., dx/2.)
		print "%.5f %.5f 2cm rlineto" % (0,-len(table12[c])/dy+dx)
		xx += (5.*dx)/3.
	print "0 0.5 moveto"
	print "0 1.5 2cm rlineto"
	print "30 0 2cm rlineto"
	print "0 -1.5 2cm rlineto"
	print "-30 0 2cm rlineto"
	print "closepath"
	print "clip"


def draw(xx, y, dx, height, last):
	if last == "*CENTROMERE*":
		print options["centromereColor"], "color"
		print "newpath"
		print "%.5f %.5f 2cm moveto" % (xx,y)
		print "%.5f %.5f 2cm rlineto" % (dx,height)
		print "%.5f %.5f 2cm rlineto" % (-dx,0)
		print "%.5f %.5f 2cm rlineto" % (dx,-height)
		print "closepath"
		print "fill"
	else:
		utils.myPsOutput.drawBox(xx, y, dx, height, last, last)

xx = 1
for c in chr1:
	
	if options["scaleY"]:
		dy = float(len(table12[c])) / 26.

	utils.myPsOutput.drawText(xx, y0, str(c), options["penColor"])
	y = y0 + 1
	
	last = ""
	nb = 0
	for (i,tmp) in table12[c]:
		if len(tmp) == 0:
			col = options["defaultColor"]
		else:
			col = tmp[0][0]
		if col == last:
			nb += 1
		else:
			if nb != 0:
				draw(xx, y, dx, nb/dy, last)
				y += nb/dy
			last = col
			nb = 1
	if nb != 0:
		draw(xx, y, dx, nb/dy, last)
	xx += (5.*dx)/3.

utils.myPsOutput.printPsFooter()
print >> sys.stderr, "OK"


