
# Module de dessin d'un karyotype

import sys
import itertools

import myTools
import myPsOutput

myTools.addModuleOptions("karyo", \
	[("dx",float,0), ("dy",float,0), ("roundedChr",bool,False), ("landscape",bool,False), ("showText",bool,True), ("drawBorder",bool,False), \
	("defaultColor",str,"black"), ("penColor",str,"black"), ("backgroundColor",str,"")] \
)

def drawKaryo(data, arguments):

	(largeur,hauteur) = myPsOutput.printPsHeader()

	if len(arguments["karyo:backgroundColor"]) > 0:
		myPsOutput.drawBox(0,0, largeur,hauteur, arguments["karyo:backgroundColor"], arguments["karyo:backgroundColor"])

	# On dessine
	dx = arguments["karyo:dx"] if arguments["karyo:dx"] > 0 else ((largeur-2.)*3.)/(5.*len(data)-2.)
	dy = arguments["karyo:dy"] if arguments["karyo:dy"] > 0 else (hauteur-4.) / float(max([len(x) for x in data.itervalues()]))
	y0 = 1.

	drawBox = myPsOutput.drawBox

	xx = 1
	for (c,l) in data.iteritems():
		def printBorder():
			print "newpath"
			print "%.5f cm %.5f cm %.5f cm 180 360 arc" % (xx+dx/2., y0+1+dx/2., dx/2.)
			print "%.5f %.5f 2cm rlineto" % (0,len(l)*dy-dx)
			print "%.5f cm %.5f cm %.5f cm 0 180 arc" % (xx+dx/2., y0+1+len(l)*dy-dx/2., dx/2.)
			print "%.5f %.5f 2cm rlineto" % (0,-len(l)*dy+dx)
			print "closepath"

		if arguments["karyo:roundedChr"]:
			print "initclip"
			printBorder()
			print "clip"

		y = y0 + 1
		for (col,items) in itertools.groupby(l):
			hauteur = len(list(items)) * dy
			if col is None:
				col = arguments["karyo:defaultColor"]
			myPsOutput.drawBox(xx, y, dx, hauteur, col, col)
			y += hauteur
		
		print "initclip"
		if arguments["karyo:drawBorder"]:
			printBorder()
			myPsOutput.setColor(arguments["karyo:penColor"], "color")
			print "stroke"

		if arguments["karyo:showText"]:
			myPsOutput.drawText(xx, y0, str(c), arguments["karyo:penColor"])
		
		xx += (5.*dx)/3.

	myPsOutput.printPsFooter()

